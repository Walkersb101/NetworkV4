#include "Quasistatic.hpp"

networkV4::quasiStaticStrain::quasiStaticStrain() {}

networkV4::quasiStaticStrain::quasiStaticStrain(
    double _maxStrain,
    StrainType _strainType,
    std::unique_ptr<BreakTypes>& _breakType)
    : m_maxStrain(_maxStrain)
    , protocol(_strainType, _breakType)
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
    bool _saveBreaks)
    : m_maxStrain(_maxStrain)
    , m_esp(_esp)
    , m_rootTol(_tol)
    , m_forceTol(_tol)
    , protocol(_strainType, _breakType)
    , m_errorOnNotSingleBreak(_errorOnNotSingleBreak)
    , m_maxStep(_maxStep)
    , m_saveBreaks(_saveBreaks)
{
}

networkV4::quasiStaticStrain::~quasiStaticStrain() {}

void networkV4::quasiStaticStrain::run(network& _network)
{
  FireMinimizer minimizer(m_forceTol);
  minimizer.integrate(_network);
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  while (true) {
    m_strainCount++;
    auto [reason, aboveThreshold] = findSingleBreak(_network);
    if (reason == SingleBreakReason::NoBreaksInStep) {
      m_dataOut->writeTimeData(genTimeData(_network, "Strain", aboveThreshold));
      m_networkOut->save(_network,
                         m_strainCount,
                         0.0,
                         "Strain-" + std::to_string(m_strainCount));
      continue;
    } else if (reason == SingleBreakReason::BreakAtLowerBound) {
      throw std::runtime_error("Break at lower bound");
    } else if (reason == SingleBreakReason::MaxStrainReached) {
      break;
    } else if (reason == SingleBreakReason::DidNotConverge) {
      throw std::runtime_error("Did not converge");
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

    if (getStrain(_network) >= m_maxStrain) {
      break;
    }
    if (m_oneBreak && breakCount > 0) {
      break;
    }
  }
}

void networkV4::quasiStaticStrain::evalStrain(network& _network,
                                              double _targetStrain,
                                              double& _maxVal,
                                              std::size_t& _count)
{
  FireMinimizer minimizer(m_forceTol);
  const double step = _targetStrain - getStrain(_network);
  strain(_network, step);
  _network.computeForces();
  minimizer.integrate(_network);
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

  network testNetwork;

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

auto networkV4::quasiStaticStrain::findSingleBreak(network& _network)
    -> std::tuple<SingleBreakReason, size_t>
{
  network testNetwork;
  double maxDistAboveA, maxDistAboveB;
  size_t breakCountA, breakCountB;
  bool converged = false;

  double a = getStrain(_network);
  double b =
      m_maxStep == 0.0 ? m_maxStrain : std::min(a + m_maxStep, m_maxStrain);

  testNetwork = _network;
  evalStrain(testNetwork, a, maxDistAboveA, breakCountA);

  if (maxDistAboveA > 0.0) {
    std::cout << "Break at lower bound with strain \n";
    std::cout << "Strain: " << a << "\n";
    std::cout << "MaxDistAboveA: " << maxDistAboveA << "\n";
    return {SingleBreakReason::BreakAtLowerBound, 0};
  }

  testNetwork = _network;
  evalStrain(testNetwork, b, maxDistAboveB, breakCountB);

  if (maxDistAboveB < 0.0) {
    _network = testNetwork;
    if (b == m_maxStrain) {
      return {SingleBreakReason::MaxStrainReached, 0};
    } else {
      return {SingleBreakReason::NoBreaksInStep, 0};
    }
  }
  try {
    converged = converge(_network,
                         a,
                         b,
                         maxDistAboveA,
                         maxDistAboveB,
                         breakCountB,
                         m_rootTol,
                         true);
  } catch (const std::invalid_argument& e) {
    evalStrain(_network, b, maxDistAboveB, breakCountB);
    return {SingleBreakReason::NoBreaksInStep, 0};
  }
  if (!converged) {
    return {SingleBreakReason::DidNotConverge, 0};
  }
  evalStrain(_network, b, maxDistAboveB, breakCountB);
  return {SingleBreakReason::Complete, breakCountB};
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
    const double force = bond.force(dist.norm());
    forces.push_back(force);
  }
  return forces;
}

auto networkV4::quasiStaticStrain::relaxBreak(network& _network) -> size_t
{
  double startEnergy, endEnergy;
  std::vector<bond*> broken;
  size_t breakCount = 0;

  OverdampedAdaptiveMinimizer integrator(m_defaultDt, m_esp);
  double nextDt = integrator.getDt();

  broken = m_breakProtocol->Break(_network, integrator);
  _network.computeForces();

  breakCount += broken.size();
  for (const auto b : broken){
    m_dataOut->writeBondData(genBondData(_network, *b));
  }

  for (size_t iter = 0; iter < m_maxMinimizerIter; ++iter) {
    startEnergy = _network.getEnergy();
    nextDt = integrator.innerIteration(_network, nextDt);
    m_t += integrator.getDt();
    endEnergy = _network.getEnergy();

    broken = m_breakProtocol->Break(_network, integrator);
    if (broken.size() > 0) {
      _network.computeForces();
      breakCount += broken.size();
      for (const auto b : broken) {
        m_dataOut->writeBondData(genBondData(_network, *b));
      }
      if (m_saveBreaks) {
        const std::string reason = "Broken-" + std::to_string(m_strainCount)
            + "-" + std::to_string(breakCount);
        m_networkOut->save(_network, m_strainCount, m_t, reason);
      }
    }

    double forceError = tools::norm(_network.getNodes().forces());
    bool converged =
        (forceError < m_forceTol) || (startEnergy - endEnergy < 1e-10);
    if (converged && broken.size() == 0) {
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
    const network& _network, const bond& _bond) -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(23 + 5 * types.size());
  data.emplace_back(m_strainCount);
  data.emplace_back(m_t);
  data.emplace_back(enumString::bondType2Str.at(_bond.type()));
  data.emplace_back(_bond.mu());
  data.emplace_back(_bond.lambda());
  data.emplace_back(_bond.naturalLength());
  data.emplace_back(_network.bondStrain(_bond));
  data.emplace_back(_bond.src());
  data.emplace_back(_bond.dst());
  data.emplace_back(_network.getNodes().position(_bond.src()).x);
  data.emplace_back(_network.getNodes().position(_bond.src()).y);
  data.emplace_back(_network.getNodes().position(_bond.dst()).x);
  data.emplace_back(_network.getNodes().position(_bond.dst()).y);
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