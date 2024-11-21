#include "Quasistatic.hpp"

networkV4::quasiStaticStrain::quasiStaticStrain() {}

networkV4::quasiStaticStrain::quasiStaticStrain(
    double _maxStrain,
    StrainType _strainType,
    std::unique_ptr<BreakTypes>& _breakType)
    : quasiStaticStrain(
          _maxStrain,
          _strainType,
          _breakType,
          config::integrators::adaptiveIntegrator::esp,
          config::rootMethods::targetTol,
          config::protocols::quasiStaticStrain::errorOnNotSingleBreak,
          _maxStrain,
          false)
{
}

networkV4::quasiStaticStrain::quasiStaticStrain(
    double _maxStrain,
    StrainType _strainType,
    std::unique_ptr<BreakTypes>& _breakType,
    double _esp,
    double _tol,
    bool _errorOnNotSingleBreak,
    double _maxStep,
    bool _single)
    : m_maxStrain(_maxStrain)
    , m_esp(_esp)
    , m_tol(_tol)
    , protocol(_strainType, _breakType)
    , m_strainCount(0)
    , m_t(0.0)
    , m_errorOnNotSingleBreak(_errorOnNotSingleBreak)
    , m_maxStep(_maxStep)
    , m_single(_single)
{
}

networkV4::quasiStaticStrain::~quasiStaticStrain() {}

void networkV4::quasiStaticStrain::run(network& _network)
{
  FireMinimizer minimizer(m_tol);
  minimizer.integrate(_network);
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  while (true) {
    m_strainCount++;
    size_t aboveThreshold = findSingleBreak(_network);
    if (aboveThreshold == 0) {
      throw std::runtime_error("No bonds above threshold");
    } else if (aboveThreshold > 1 && m_errorOnNotSingleBreak) {
      throw std::runtime_error("More than one bond above threshold");
    }
    m_t = 0.0;

    m_dataOut->writeTimeData(genTimeData(_network, "Start", aboveThreshold));
    m_networkOut->save(
        _network, m_strainCount, 0.0, "Start-" + std::to_string(m_strainCount));

    size_t breakCount = relaxBreak(_network);

    m_dataOut->writeTimeData(genTimeData(_network, "End", breakCount));
    m_networkOut->save(
        _network, m_strainCount, m_t, "End-" + std::to_string(m_strainCount));

    if (m_single || getStrain(_network) >= m_maxStrain) {
        break;
    }
  }
}

void networkV4::quasiStaticStrain::evalStrain(network& _network,
                                              double _step,
                                              double& _maxVal,
                                              std::size_t& _count,
                                              bool _save)
{
  FireMinimizer minimizer(config::integrators::miminizer::tol);

  double targetStrain = getStrain(_network) + _step;
  while (getStrain(_network) < targetStrain - 1e-15) {
    const double strainStep =
        std::min(targetStrain - getStrain(_network), m_maxStep);
    strain(_network, strainStep);
    _network.computeForces();
    minimizer.integrate(_network);
    if (_save) {
      m_dataOut->writeTimeData(genTimeData(_network, "Strain", 0));
    }
  }
  m_breakProtocol->Data(_network, minimizer, _maxVal, _count);
}

auto networkV4::quasiStaticStrain::converge(network& _baseNetwork,
                                            double& _a,
                                            double& _b,
                                            double& _maxDistAboveA,
                                            double& _maxDistAboveB,
                                            size_t& _breakCountB,
                                            double _tol,
                                            bool _earlyExit) -> bool
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
    } else {
      _a = xITP;
      _maxDistAboveA = maxDistAboveITP;
    }

    if (std::abs(_b - _a) < 2 * _tol || (_earlyExit && _breakCountB == 1)) {
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
  network testNetwork;

  testNetwork = _network;
  evalStrain(testNetwork, a, maxDistAboveA, breakCountA);

  if (maxDistAboveA >= 0.0) {
    return 0;
  }

  double guessStrain =
      std::min(std::abs(maxDistAboveA)
                   * config::protocols::quasiStaticStrain::strainGuessScale,
               b);

  testNetwork = _network;
  evalStrain(testNetwork, guessStrain, maxDistAboveB, breakCountB);

  if (maxDistAboveA * maxDistAboveB > 0.0) {  // guess was not big enough
    a = guessStrain;
    maxDistAboveA = maxDistAboveB;
    testNetwork = _network;
    evalStrain(testNetwork, b, maxDistAboveB, breakCountB);
  } else {
    b = guessStrain;
  }

  bool converged = converge(
      _network, a, b, maxDistAboveA, maxDistAboveB, breakCountB, m_tol);
  if (converged && breakCountB != 1) {
    converged = converge(_network,
                         a,
                         b,
                         maxDistAboveA,
                         maxDistAboveB,
                         breakCountB,
                         config::rootMethods::minTol,
                         true);
  }
  evalStrain(_network, b, maxDistAboveB, breakCountB, true);
  return breakCountB;
}

std::vector<double> networkV4::forceMags(const network& _network)
{
  const auto& bonds = _network.getBonds();
  const auto& position = _network.getNodes().positions();

  std::vector<double> forces;
  forces.reserve(bonds.size());
  for (const auto& bond : bonds) {
    const vec2d dist =
        _network.minDist(position[bond.src()], position[bond.dst()]);
    const double force = bond.force(dist.length());
    forces.push_back(force);
  }
  return forces;
}

auto networkV4::quasiStaticStrain::relaxBreak(network& _network) -> size_t
{
  std::size_t maxIter = config::integrators::miminizer::maxIter;

  OverdampedAdaptiveMinimizer integrator(config::integrators::default_dt,
                                         m_esp);
  double nextDt = integrator.getDt();

  std::vector<size_t> broken = m_breakProtocol->Break(_network, integrator);
  size_t breakCount = broken.size();
  if (breakCount == 0) {
    return 0;
  }
  for (const auto& b : broken) {
    m_dataOut->writeBondData(genBondData(_network, b));
  }

  for (size_t iter = 0; iter < maxIter; ++iter) {
    nextDt = integrator.innerIteration(_network, nextDt);
    m_t += integrator.getDt();

    broken = m_breakProtocol->Break(_network, integrator);
    breakCount += broken.size();

    for (const auto& b : broken) {
      m_dataOut->writeBondData(genBondData(_network, b));
    }

    if (broken.size() > 0) {
      m_networkOut->save(_network,
                         m_strainCount,
                         m_t,
                         "Broken-" + std::to_string(m_strainCount) + "-"
                             + std::to_string(breakCount));
    }

    double error = tools::norm(_network.getNodes().forces());
    if (error < m_tol && broken.size() == 0) {
      return breakCount;
    }
    // std::cout << std::setprecision(10) << iter << " " << error << " "
    //           << _network.getEnergy() << " " << integrator.getDt() << " "
    //           << tools::maxLength(_network.getNodes().forces()) << " "
    //           << tools::maxAbsComponent(_network.getNodes().forces())
    //           << std::endl;
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
      enumString::strainType2Str.at(m_strainType) + "Strain",
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
      enumString::strainType2Str.at(m_strainType) + "Strain",
      "RMSForce",
      "MaxForce"};

  tensorData::addBondCountByTypeHeader(dataHeader, types);
  tensorData::addBondCountByTypeHeader(bondHeader, types);

  tensorData::addStressByTypeHeader(bondHeader, types);
  tensorData::addStressByTypeHeader(dataHeader, types);

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
  data.emplace_back(_reason);
  data.emplace_back(m_strainCount);
  data.emplace_back(_breakCount);
  data.emplace_back(m_t);
  data.emplace_back(_network.getDomain().x);
  data.emplace_back(_network.getDomain().y);
  data.emplace_back(getStrain(_network));

  data.emplace_back(_network.getBonds().connectedCount());
  for (const auto& type : types) {
    data.emplace_back(_network.getBonds().connectedCount(type));
  }

  data.insert(
      data.end(),
      {globalStress.xx, globalStress.xy, globalStress.yx, globalStress.yy});
  for (const auto& type : types) {
    tensor2 stress = stresses[type];
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
  data.emplace_back(m_strainCount);
  data.emplace_back(m_t);
  data.emplace_back(enumString::bondType2Str.at(b.type()));
  data.emplace_back(b.mu());
  data.emplace_back(b.lambda());
  data.emplace_back(b.naturalLength());
  data.emplace_back(_network.bondStrain(b));
  data.emplace_back(b.src());
  data.emplace_back(b.dst());
  data.emplace_back(_network.getNodes().position(b.src()).x);
  data.emplace_back(_network.getNodes().position(b.src()).y);
  data.emplace_back(_network.getNodes().position(b.dst()).x);
  data.emplace_back(_network.getNodes().position(b.dst()).y);
  data.emplace_back(_network.getDomain().x);
  data.emplace_back(_network.getDomain().y);
  data.emplace_back(getStrain(_network));
  data.emplace_back(tools::norm(_network.getNodes().forces()));
  data.emplace_back(tools::maxLength(_network.getNodes().forces()));
  data.emplace_back(_network.getBonds().connectedCount());
  for (const auto& type : types) {
    data.emplace_back(_network.getBonds().connectedCount(type));
  }

  data.insert(
      data.end(),
      {globalStress.xx, globalStress.xy, globalStress.yx, globalStress.yy});
  for (const auto& type : types) {
    tensor2 stress = stresses[type];
    data.insert(data.end(), {stress.xx, stress.xy, stress.yx, stress.yy});
  }
  return data;
}