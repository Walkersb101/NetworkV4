#include <algorithm>
#include <cassert>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Protocol.hpp"

#include "DataOut.hpp"
#include "Intergrator.hpp"
#include "Network.hpp"
#include "Roots.hpp"
#include "Tools.hpp"

networkV4::protocol::protocol() {}

networkV4::protocol::protocol(StrainType _strainType, BreakType _breakType)
    : m_strainType(_strainType)
    , m_breakType(_breakType)
    , m_dataOut(nullptr)
{
}

networkV4::protocol::~protocol() {}

double networkV4::protocol::getStrain(const network& _network) const
{
  switch (m_strainType) {
    case StrainType::Shear:
      return _network.getShearStrain();
    case StrainType::Elongation:
      return _network.getElongationStrain();
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::protocol::strain(network& _network, double _step)
{
  switch (m_strainType) {
    case StrainType::Shear:
      _network.shear(_step);
      break;
    case StrainType::Elongation:
      _network.elongate(_step);
      break;
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::protocol::breakData(const network& _network,
                                    BreakType _type,
                                    double& _maxDistAbove,
                                    std::size_t& _count)
{
  switch (_type) {
    case BreakType::None:
      _maxDistAbove = 0.0;
      _count = 0;
      break;
    case BreakType::Strain:
      _network.strainBreakData(_maxDistAbove, _count);
      break;
    default:
      throw std::runtime_error("Unknown break type");
  }
}

auto networkV4::protocol::breakBonds(network& _network, BreakType _type) const
    -> std::vector<std::size_t>
{
  switch (_type) {
    case BreakType::None:
      return std::vector<std::size_t>();
    case BreakType::Strain:
      return _network.strainBreak();
    default:
      throw std::runtime_error("Unknown break type");
  }
}

void networkV4::protocol::run(network& _network) {}

void networkV4::protocol::initIO(const network& _network,
                                 std::unique_ptr<dataOut>& _dataOut,
                                 std::unique_ptr<networkOut>& _networkOut)
{
}

networkV4::quasiStaticStrain::quasiStaticStrain() {}

networkV4::quasiStaticStrain::quasiStaticStrain(double _maxStrain,
                                                StrainType _strainType,
                                                BreakType _breakType)
    : m_maxStrain(_maxStrain)
    , protocol(_strainType, _breakType)
    , m_esp(config::intergrators::adaptiveIntergrator::esp)
    , m_tol(config::rootMethods::targetTol)
    , m_strainCount(0)
    , m_t(0.0)
    , m_errorOnNotSingleBreak(
          config::protocols::quasiStaticStrain::errorOnNotSingleBreak)
{
}

networkV4::quasiStaticStrain::quasiStaticStrain(double _maxStrain,
                                                StrainType _strainType,
                                                BreakType _breakType,
                                                double _esp,
                                                double _tol,
                                                bool _errorOnNotSingleBreak)
    : m_maxStrain(_maxStrain)
    , m_esp(_esp)
    , m_tol(_tol)
    , protocol(_strainType, _breakType)
    , m_strainCount(0)
    , m_t(0.0)
    , m_errorOnNotSingleBreak(_errorOnNotSingleBreak)
{
}

networkV4::quasiStaticStrain::~quasiStaticStrain() {}

void networkV4::quasiStaticStrain::run(network& _network)
{
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  while (true) {
    size_t aboveThreshold = findSingleBreak(_network);
    if (aboveThreshold == 0) {
      std::runtime_error("No bonds above threshold");
      return;
    } else if (aboveThreshold > 1 && m_errorOnNotSingleBreak) {
      std::runtime_error("More than one bond above threshold");
      return;
    }
    m_strainCount++;
    m_t = 0.0;

    m_dataOut->writeTimeData(genTimeData(_network, "Start", aboveThreshold));
    m_networkOut->save(_network, "Start-" + std::to_string(m_strainCount));

    size_t breakCount = relaxBreak(_network);

    m_dataOut->writeTimeData(genTimeData(_network, "End", breakCount));
    m_networkOut->save(_network, "End-" + std::to_string(m_strainCount));
  }
}

void networkV4::quasiStaticStrain::evalStrain(network& _network,
                                              double _step,
                                              double& _maxVal,
                                              std::size_t& _count)
{
  // FireMinimizer minimizer(config::miminizer::tol);
  OverdampedAdaptiveMinimizer minimizer(
      config::intergrators::default_dt,
      config::intergrators::adaptiveIntergrator::esp,
      config::intergrators::miminizer::tol);
  strain(_network, _step);
  minimizer.integrate(_network);
  breakData(_network, m_breakType, _maxVal, _count);
}

auto networkV4::quasiStaticStrain::converge(network& _baseNetwork,
                                            network& _networkB,
                                            double& _a,
                                            double& _b,
                                            double& _maxDistAboveA,
                                            double& _maxDistAboveB,
                                            size_t& _breakCountB,
                                            double _tol) -> bool
{
  if (_maxDistAboveA * _maxDistAboveB > 0.0 || _a >= _b) {
    std::invalid_argument("Root not bracketed");
    return false;
  }
  roots::ITP solver(_a, _b, _tol);

  network testNetwork = _baseNetwork;

  for (std::size_t iters = 0; iters < solver.nMax(); ++iters) {
    double xITP = solver.guessRoot(_a, _b, _maxDistAboveA, _maxDistAboveB);

    double maxDistAboveITP;
    std::size_t breakCountITP;
    testNetwork = _baseNetwork;
    evalStrain(testNetwork, xITP, maxDistAboveITP, breakCountITP);
    if (maxDistAboveITP >= 0.) {
      _b = xITP;
      _maxDistAboveB = maxDistAboveITP;
      _breakCountB = breakCountITP;
      _networkB = testNetwork;
    } else {
      _a = xITP;
      _maxDistAboveA = maxDistAboveITP;
    }

    if (std::abs(_b - _a) < 2 * m_tol) {
      return true;
    }
  }
  return false;
}

auto networkV4::quasiStaticStrain::findSingleBreak(network& _network) -> size_t
{
  double a = 0.0;
  double b = m_maxStrain - getStrain(_network);
  if (b < 0.0) {
    return 0;
  }

  double maxDistAboveA, maxDistAboveB;
  std::size_t breakCountA, breakCountB;

  network testNetwork = _network;
  evalStrain(testNetwork, a, maxDistAboveA, breakCountA);
  network networkB = _network;
  evalStrain(networkB, b, maxDistAboveB, breakCountB);

  // TODO: Print the results
  bool converged = false;

  converged = converge(_network,
                       networkB,
                       a,
                       b,
                       maxDistAboveA,
                       maxDistAboveB,
                       breakCountB,
                       m_tol);
  if (converged && breakCountB != 1) {
    converged = converge(_network,
                         networkB,
                         a,
                         b,
                         maxDistAboveA,
                         maxDistAboveB,
                         breakCountB,
                         config::rootMethods::minTol);
  }
  _network = networkB;
  return converged ? breakCountB : 0;
}

auto networkV4::quasiStaticStrain::relaxBreak(network& _network) -> size_t
{
  std::size_t maxIter = config::intergrators::miminizer::maxIter;

  OverdampedAdaptiveEulerHeun integrator;

  std::vector<size_t> broken = breakBonds(_network, m_breakType);
  size_t breakCount = broken.size();
  if (breakCount == 0) {
    return 0;
  }
  for (const auto& b : broken) {
    m_dataOut->writeBondData(genBondData(_network, b));
  }

  for (size_t iter = 0; iter < maxIter; ++iter) {
    integrator.integrate(_network);

    broken = breakBonds(_network, m_breakType);
    breakCount += broken.size();

    for (const auto& b : broken) {
      m_dataOut->writeBondData(genBondData(_network, b));
    }

    double error = tools::maxAbsComponent(_network.getNodes().forces());
    if (error < m_tol && broken.size() == 0) {
      break;
    }

    m_t += integrator.getDt();
  }
  return breakCount;
}

void networkV4::quasiStaticStrain::initIO(
    const network& _network,
    std::unique_ptr<dataOut>& _dataOut,
    std::unique_ptr<networkOut>& _networkOut)
{
  auto types = tools::uniqueBondTypes(_network.getBonds());

  std::vector<std::string> dataHeader = {
      "Reason",
      "StrainCount",
      "BreakCount",
      "Time",
      "Domainx",
      "Domainy",
      enum2str::strainTypeString(m_strainType) + "Strain",
  };

  std::vector<std::string> bondHeader = {
      "StrainCount",
      "Time",
      "Type",
      "mu",
      "lambda",
      "NaturalLength",
      "BondStrain",
      "Index1",
      "Index2",
      "x1",
      "y1",
      "x2",
      "y2",
      "Domainx",
      "Domainy",
      enum2str::strainTypeString(m_strainType) + "Strain",
      "RMSForce",
      "MaxForce"};

  enum2str::addBondCountByTypeHeader(dataHeader, types);
  enum2str::addBondCountByTypeHeader(bondHeader, types);

  enum2str::addStressByTypeHeader(bondHeader, types);
  enum2str::addStressByTypeHeader(dataHeader, types);

  _dataOut->initFiles(dataHeader, bondHeader);
  m_dataOut = _dataOut.get();
  m_networkOut = _networkOut.get();
}

auto networkV4::quasiStaticStrain::genTimeData(const network& _network,
                                               const std::string& _reason,
                                               size_t _breakCount)
    -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(12 + 5 * types.size());
  data.push_back(_reason);
  data.push_back(m_strainCount);
  data.push_back(_breakCount);
  data.push_back(m_t);
  data.push_back(_network.getDomain().x);
  data.push_back(_network.getDomain().y);
  data.push_back(getStrain(_network));

  data.push_back(_network.getBonds().connectedCount());
  for (const auto& type : types) {
    data.push_back(_network.getBonds().connectedCount(type));
  }

  data.insert(
      data.end(),
      {globalStress.xx, globalStress.xy, globalStress.yx, globalStress.yy});
  for (const auto& type : types) {
    tensor2 stress = stresses.at(type);
    data.insert(data.end(), {stress.xx, stress.xy, stress.yx, stress.yy});
  }
  return data;
}

auto networkV4::quasiStaticStrain::genBondData(
    const network& _network, size_t _bondIndex) -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();
  auto& b = _network.getBonds().get(_bondIndex);

  std::vector<writeableTypes> data;
  data.reserve(23 + 5 * types.size());
  data.push_back(m_strainCount);
  data.push_back(m_t);
  data.push_back(enum2str::bondTypeString(b.type()));
  data.push_back(b.mu());
  data.push_back(b.lambda());
  data.push_back(b.naturalLength());
  data.push_back(_network.bondStrain(b));
  data.push_back(b.src());
  data.push_back(b.dst());
  data.push_back(_network.getNodes().position(b.src()).x);
  data.push_back(_network.getNodes().position(b.src()).y);
  data.push_back(_network.getNodes().position(b.dst()).x);
  data.push_back(_network.getNodes().position(b.dst()).y);
  data.push_back(_network.getDomain().x);
  data.push_back(_network.getDomain().y);
  data.push_back(getStrain(_network));
  data.push_back(tools::norm(_network.getNodes().forces()));
  data.push_back(tools::maxLength(_network.getNodes().forces()));
  data.push_back(_network.getBonds().connectedCount());
  for (const auto& type : types) {
    data.push_back(_network.getBonds().connectedCount(type));
  }

  data.insert(
      data.end(),
      {globalStress.xx, globalStress.xy, globalStress.yx, globalStress.yy});
  for (const auto& type : types) {
    tensor2 stress = stresses.at(type);
    data.insert(data.end(), {stress.xx, stress.xy, stress.yx, stress.yy});
  }
  return data;
}