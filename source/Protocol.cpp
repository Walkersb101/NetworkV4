#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Protocol.hpp"

#include "DataOut.hpp"
#include "Intergrator.hpp"
#include "Network.hpp"
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
    , m_esp(config::adaptiveIntergrator::esp)
    , m_tol(config::ITPMethod::tol)
    , m_strainCount(0)
    , m_t(0.0)
{
}

networkV4::quasiStaticStrain::quasiStaticStrain(double _maxStrain,
                                                StrainType _strainType,
                                                BreakType _breakType,
                                                double _esp,
                                                double _tol)
    : m_maxStrain(_maxStrain)
    , m_esp(_esp)
    , m_tol(_tol)
    , protocol(_strainType, _breakType)
    , m_strainCount(0)
    , m_t(0.0)
{
}

networkV4::quasiStaticStrain::~quasiStaticStrain() {}

void networkV4::quasiStaticStrain::run(network& _network)
{
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  while (true) {
    if (!findSingleBreak(_network)) {
      return;
    }
    m_strainCount++;
    m_t = 0.0;
    m_networkOut->save(_network, "Start-" + std::to_string(m_strainCount)); 
    m_dataOut->writeTimeData(genTimeData(_network, "Start", 0));
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
  OverdampedAdaptiveMinimizer minimizer(config::default_dt,
                                        config::adaptiveIntergrator::esp,
                                        config::miminizer::tol);
  strain(_network, _step);
  minimizer.integrate(_network);
  breakData(_network, m_breakType, _maxVal, _count);
}

double sign(const double _val)
{
  return _val > 0.0 ? 1.0 : -1.0;
}

auto networkV4::quasiStaticStrain::findSingleBreak(network& _network) -> bool
{
  size_t n0 = config::ITPMethod::n0;
  double k1Scale = config::ITPMethod::k1Scale;
  double k2 = config::ITPMethod::k2;

  double strain = getStrain(_network);
  double a = 0.0;
  double b = m_maxStrain - strain;

  if (b < 0.0) {
    return false;
  }

  double maxDistAboveA, maxDistAboveB;
  std::size_t breakCountA, breakCountB;

  network testNetwork = _network;
  evalStrain(testNetwork, 0.0, maxDistAboveA, breakCountA);
  network networkB = _network;
  evalStrain(networkB, b, maxDistAboveB, breakCountB);

  // TODO: Print the results

  if (maxDistAboveA * maxDistAboveB > 0.0) {
    // TODO: print Error
    return false;
  }

  size_t nhalf = std::ceil(std::log2((b - a) / (2 * m_tol)));
  size_t nmax = nhalf + n0;

  double k1 = k1Scale * (b - a);

  for (std::size_t iters = 0; iters < nmax; ++iters) {
    double xhalf = (a + b) * 0.5;
    double r = m_esp * std::pow(2, nmax - iters) - (b - a) * 0.5;
    double delta = k1 * std::pow(b - a, k2);

    double alpha = iters % 3 == 0 ? 0.5 : 1.0;
    double xf = (alpha * maxDistAboveB * a - maxDistAboveA * b)
        / (alpha * maxDistAboveB - maxDistAboveA);
    double sigma = sign(xhalf - xf);
    double xt = delta < std::abs(xf - xhalf) ? xhalf + sigma * delta : xhalf;
    double xITP = std::abs(xt - xhalf) <= r ? xt : xhalf + sigma * r;

    double maxDistAboveITP;
    std::size_t breakCountITP;
    testNetwork = _network;
    evalStrain(testNetwork, xITP, maxDistAboveITP, breakCountITP);
    if (maxDistAboveITP > 0.) {
      b = xITP;
      maxDistAboveB = maxDistAboveITP;
      breakCountB = breakCountITP;
      networkB = testNetwork;
    } else if (maxDistAboveITP < 0.) {
      a = xITP;
      maxDistAboveA = maxDistAboveITP;
      breakCountA = breakCountITP;
    } else {
      _network = testNetwork;
      return true;
    }

    if (std::abs(b - a) < 2 * m_tol && breakCountB == 1) {
      _network = networkB;
      return true;
    }
  }

  return false;
}

auto networkV4::quasiStaticStrain::relaxBreak(network& _network) -> size_t
{
  std::size_t maxIter = config::miminizer::maxIter;

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