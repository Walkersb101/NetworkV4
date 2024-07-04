#include <algorithm>
#include <cmath>
#include <numeric>
#include <stdexcept>

#include "Protocol.hpp"

#include "DataOut.hpp"
#include "Intergrator.hpp"
#include "Network.hpp"

auto uniqueBondTypes(const networkV4::bonds& _bonds)
    -> std::vector<networkV4::bondType>
{
  std::vector<networkV4::bondType> types;
  for (const auto& bond : _bonds) {
    if (std::find(types.begin(), types.end(), bond.type()) == types.end()) {
      types.push_back(bond.type());
    }
  }
  return types;
}

auto typeString(const networkV4::bondType& _type) -> std::string
{
  switch (_type) {
    case networkV4::bondType::single:
      return "single";
    case networkV4::bondType::matrix:
      return "matrix";
    case networkV4::bondType::sacrificial:
      return "sacrificial";
    default:
      throw std::runtime_error("Unknown bond type");
  }
}

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
      break;
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
      break;
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
      break;
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
      break;
  }
}

void networkV4::protocol::run(network& _network) {}

void networkV4::protocol::initIO(const network& _network,
                                 std::unique_ptr<dataOut>& _dataOut)
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
{
}

networkV4::quasiStaticStrain::~quasiStaticStrain() {}

void networkV4::quasiStaticStrain::run(network& _network)
{
  while (true) {
    bool findBreak = findSingleBreak(_network);
    if (!findBreak) {
      return;
    }
    m_t = 0.0;
    size_t breakCount = relaxBreak(_network);
  }
}

void networkV4::quasiStaticStrain::evalStrain(network& _network,
                                              double _step,
                                              double& _maxVal,
                                              std::size_t& _count)
{
  FireMinimizer minimizer(config::miminizer::tol);
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
  std::size_t iters = 0;

  while (iters < nmax) {
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
      a = xITP;
      b = xITP;
      maxDistAboveA = 0.0;
      maxDistAboveB = 0.0;
      breakCountA = breakCountITP;
      breakCountB = breakCountITP;
      networkB = testNetwork;
    }

    if (std::abs(b - a) < 2 * m_tol && breakCountB == 1) {
      _network = networkB;
      return true;
    }
    iters++;
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
  std::size_t iters = 0;
  while (iters < maxIter) {
    integrator.integrate(_network);

    broken = breakBonds(_network, m_breakType);
    breakCount += broken.size();

    double error = std::accumulate(_network.getNodes().forces().begin(),
                                   _network.getNodes().forces().end(),
                                   0.0,
                                   [](double _max, const vec2d& _v)
                                   { return std::max(_max, _v.abs().max()); });
    if (error < m_tol && broken.size() == 0) {
      break;
    }

    m_t += integrator.getDt();
    ++iters;
  }
  return breakCount;
}


void networkV4::quasiStaticStrain::initIO(const network& _network,
                                          std::unique_ptr<dataOut>& _dataOut)
{
  auto types = uniqueBondTypes(_network.getBonds());

  std::vector<std::string> dataHeader = {"Reason",
                                         "Time",
                                         "dt",
                                         "Domainx",
                                         "Domainy",
                                         "ElongationStrain",
                                         "ShearStrain",
                                         "RMSVel",
                                         "MaxVel"};

  std::vector<std::string> bondHeader = {"Time",
                                         "dt",
                                         "BreakCount",
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
                                         "ElongationStrain",
                                         "ShearStrain",
                                         "RMSVel",
                                         "MaxVel"};

  dataHeader.push_back("BondCount");
  bondHeader.push_back("BondCount");
  for (const auto& type : types) {
    dataHeader.push_back(typeString(type) + "Count");
    bondHeader.push_back(typeString(type) + "Count");
  }
  dataHeader.insert(dataHeader.end(),
                    {"Stressxx", "Stressxy", "Stressyx", "Stressyy"});
  bondHeader.insert(bondHeader.end(),
                    {"Stressxx", "Stressxy", "Stressyx", "Stressyy"});
  for (const auto& type : types) {
    dataHeader.insert(dataHeader.end(),
                      {typeString(type) + "Stressxx",
                       typeString(type) + "Stressxy",
                       typeString(type) + "Stressyx",
                       typeString(type) + "Stressyy"});
    bondHeader.insert(bondHeader.end(),
                      {typeString(type) + "Stressxx",
                       typeString(type) + "Stressxy",
                       typeString(type) + "Stressyx",
                       typeString(type) + "Stressyy"});
  }

  _dataOut->initFiles(dataHeader, bondHeader);
  m_dataOut = _dataOut.get();
}
