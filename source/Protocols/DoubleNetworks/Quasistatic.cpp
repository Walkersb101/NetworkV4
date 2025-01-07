#include <iostream>

#include "Quasistatic.hpp"

networkV4::protocols::quasiStaticStrainDouble::quasiStaticStrainDouble(
    std::shared_ptr<deform::deformBase>& _deform,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
    std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
    const network& _network,
    double _maxStrain,
    double _rootTol,
    integration::AdaptiveParams _params,
    const minimisation::minimiserParams& _minParams,
    bool _errorOnNotSingleBreak,
    double _maxStep,
    bool _saveBreaks)
    : protocolBase(_deform, _dataOut, _bondsOut, _networkOut, _network)
    , m_maxStrain(_maxStrain)
    , m_rootTol(_rootTol)
    , m_params(_params)
    , m_minParams(_minParams)
    , m_errorOnNotSingleBreak(_errorOnNotSingleBreak)
    , m_maxStep(_maxStep)
    , m_saveBreaks(_saveBreaks)
{
  std::vector<IO::timeSeries::writeableTypes> dataHeader = {
      "Reason",
      "StrainCount",
      "BreakCount",
      "Time",
      "Domainx",
      "Domainy",
      "ShearStrain",
      "ElongationStrainX",
      "ElongationStrainY",
      "BondCount",
      "SacrificialBondCount",
      "MatrixBondCount",
      "StressXX",
      "StressXY",
      "StressYX",
      "StressYY",
      "MatrixStressXX",
      "MatrixStressXY",
      "MatrixStressYX",
      "MatrixStressYY",
      "SacrificialStressXX",
      "SacrificialStressXY",
      "SacrificialStressYX",
      "SacrificialStressYY"};

  std::vector<IO::timeSeries::writeableTypes> bondHeader = {
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
      "ElongationStrainX",
      "ElongationStrainY",
      //"RMSForce",
      //"MaxForce"
      "BondCount",
      "SacrificialBondCount",
      "MatrixBondCount",
      "StressXX",
      "StressXY",
      "StressYX",
      "StressYY",
      "MatrixStressXX",
      "MatrixStressXY",
      "MatrixStressYX",
      "MatrixStressYY",
      "SacrificialStressXX",
      "SacrificialStressXY",
      "SacrificialStressYX",
      "SacrificialStressYY"};

  _dataOut->write(dataHeader);
  _bondsOut->write(bondHeader);
}

networkV4::protocols::quasiStaticStrainDouble::~quasiStaticStrainDouble() {}

void networkV4::protocols::quasiStaticStrainDouble::run(network& _network)
{
  minimisation::fire2 minimizer(m_minParams);

  minimizer.minimise(_network);
  m_dataOut->write(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  while (true) {
    m_strainCount++;
    auto [reason, aboveThreshold] = findSingleBreak(_network);
    if (reason == SingleBreakReason::NoBreaksInStep) {
      m_dataOut->write(genTimeData(_network, "Strain", aboveThreshold));
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

    m_dataOut->write(genTimeData(_network, "Start", aboveThreshold));
    m_networkOut->save(
        _network, m_strainCount, 0.0, "Start-" + std::to_string(m_strainCount));

    size_t breakCount = relaxBreak(_network);

    m_dataOut->write(genTimeData(_network, "End", breakCount));
    m_networkOut->save(
        _network, m_strainCount, m_t, "End-" + std::to_string(m_strainCount));

    if (m_deform->getStrain(_network) >= m_maxStrain) {
      break;
    }
    if (m_oneBreak && breakCount > 0) {
      break;
    }
  }
}

auto networkV4::protocols::quasiStaticStrainDouble::evalStrain(
    network& _network, double _targetStrain) -> std::tuple<double, size_t>
{
  minimisation::fire2 minimizer(m_minParams);
  const double step = _targetStrain - m_deform->getStrain(_network);
  m_deform->strain(_network, step);
  minimizer.minimise(_network);
  return breakData(_network);
}

auto networkV4::protocols::quasiStaticStrainDouble::converge(
    network& _baseNetwork,
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

    testNetwork = _baseNetwork;
    auto [maxDistAboveITP, breakCountITP] = evalStrain(testNetwork, xITP);
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

auto networkV4::protocols::quasiStaticStrainDouble::findSingleBreak(
    network& _network) -> std::tuple<SingleBreakReason, size_t>
{
  network testNetwork = _network;
  double maxDistAboveA, maxDistAboveB;
  size_t breakCountA, breakCountB;
  bool converged = false;

  double a = m_deform->getStrain(_network);
  double b = std::min(a + m_maxStep, m_maxStrain);

  testNetwork = _network;
  std::tie(maxDistAboveA, breakCountA) = evalStrain(testNetwork, a);

  // TODO : change to Log
  if (maxDistAboveA > 0.0) {
    std::cout << "Break at lower bound with strain \n";
    std::cout << "Strain: " << a << "\n";
    std::cout << "MaxDistAboveA: " << maxDistAboveA << "\n";
    return {SingleBreakReason::BreakAtLowerBound, 0};
  }

  testNetwork = _network;
  std::tie(maxDistAboveB, breakCountB) = evalStrain(testNetwork, b);

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
    std::tie(maxDistAboveB, breakCountB) = evalStrain(_network, b);
    return {SingleBreakReason::NoBreaksInStep, 0};
  }
  if (!converged) {
    return {SingleBreakReason::DidNotConverge, 0};
  }
  std::tie(maxDistAboveB, breakCountB) = evalStrain(_network, b);
  return {SingleBreakReason::Complete, breakCountB};
}

auto networkV4::protocols::quasiStaticStrainDouble::relaxBreak(
    network& _network) -> size_t
{
  double Ecurr = _network.getEnergy();
  double Eprev = Ecurr;
  size_t breakCount = 0;
  m_t = 0.0;

  // TODO: allow changing Gamma
  integration::AdaptiveOverdampedEulerHeun stepper(1.0, m_params);

  _network.computeForces(true);

  breakCount += _network.getBreakQueue().size();
  while (!_network.getBreakQueue().empty()) {
    const auto& broken = _network.getBreakQueue().front();
    m_bondsOut->write(genBondData(_network, broken));
    _network.getBreakQueue().pop_front();
  }

  size_t iter = 0;
  while (iter++ < m_minParams.maxIter) {
    Eprev = Ecurr;
    stepper.step(_network);
    _network.computeForces(true);
    Ecurr = _network.getEnergy();
    m_t += stepper.getDt();

    bool brokenInStep = !_network.getBreakQueue().empty();
    breakCount += _network.getBreakQueue().size();
    if (brokenInStep && m_saveBreaks) {
      const std::string reason = "Broken-" + std::to_string(m_strainCount) + "-"
          + std::to_string(breakCount);
      m_networkOut->save(_network, m_strainCount, m_t, reason);
    }

    while (!_network.getBreakQueue().empty()) {
      const auto& broken = _network.getBreakQueue().front();
      m_bondsOut->write(genBondData(_network, broken));
      _network.getBreakQueue().pop_front();
    }

    if (!brokenInStep
        && fabs(Ecurr - Eprev) < m_minParams.Etol * 0.5
                * (fabs(Ecurr) + fabs(Eprev) + minimisation::EPS_ENERGY))
    {
      return breakCount;
    }

    double fdotf =
        xdoty(_network.getNodes().forces(), _network.getNodes().forces());
    if (!brokenInStep && fdotf < m_minParams.Ftol * m_minParams.Ftol) {
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

auto networkV4::protocols::quasiStaticStrainDouble::genTimeData(
    const network& _network,
    const std::string& _reason,
    size_t _breakCount) -> std::vector<IO::timeSeries::writeableTypes>
{
  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress = stresses.get(_network.getTags().get("sacrificial"));
  const auto& matStress = stresses.get(_network.getTags().get("matrix"));

  const auto& box = _network.getBox();
  auto counts = getCounts(_network);

  return {_reason,
          m_strainCount,
          _breakCount,
          m_t,
          box.getLx(),
          box.getLy(),
          box.shearStrain(),
          _network.getElongationStrain().x,
          _network.getElongationStrain().y,
          std::get<0>(counts),
          std::get<1>(counts),
          std::get<2>(counts),
          globalStress.xx,
          globalStress.xy,
          globalStress.yx,
          globalStress.yy,
          matStress.xx,
          matStress.xy,
          matStress.yx,
          matStress.yy,
          sacStress.xx,
          sacStress.xy,
          sacStress.yx,
          sacStress.yy};
}

auto networkV4::protocols::quasiStaticStrainDouble::genBondData(
    const network& _network,
    const breakInfo& _bond) -> std::vector<IO::timeSeries::writeableTypes>
{
  auto sacTag = _network.getTags().get("sacrificial");
  auto matTag = _network.getTags().get("matrix");

  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress = stresses.get(sacTag);
  const auto& matStress = stresses.get(matTag);

  const auto& box = _network.getBox();
  auto counts = getCounts(_network);

  const auto& nodes = _network.getNodes();
  const auto bonds = _network.getBonds();

  const auto& binfo = std::get<0>(_bond);
  const auto& type = std::get<1>(_bond);
  const auto& brk = std::get<2>(_bond);

  bool harmonic = std::holds_alternative<Forces::HarmonicBond>(type);
  bool strainBreak = std::holds_alternative<BreakTypes::StrainBreak>(brk);
  bool sacrificial = _network.getBondTags().hasTag(binfo.index, sacTag);

  const auto& pos1 = nodes.positions()[binfo.src];
  const auto& pos2 = nodes.positions()[binfo.dst];
  size_t bondSrc = nodes.indices()[binfo.src];
  size_t bondDst = nodes.indices()[binfo.dst];

  double r = box.minDist(pos1, pos2).norm();
  double r0 = harmonic ? std::get<Forces::HarmonicBond>(type).r0() : 0.0;

  return {
      m_strainCount,
      m_t,
      sacrificial ? "Sacrificial" : "Matrix",
      harmonic ? std::get<Forces::HarmonicBond>(std::get<1>(_bond)).k() : 0.0,
      strainBreak ? std::get<BreakTypes::StrainBreak>(brk).lambda() : 0.0,
      harmonic ? std::get<Forces::HarmonicBond>(std::get<1>(_bond)).r0() : 0.0,
      harmonic ? (r - r0) / r0 : 0.0,
      bondSrc,
      bondDst,
      pos1.x,
      pos1.y,
      pos2.x,
      pos2.y,
      box.getLx(),
      box.getLy(),
      _network.getElongationStrain().x,
      _network.getElongationStrain().y,
      std::get<0>(counts),
      std::get<1>(counts),
      std::get<2>(counts),
      globalStress.xx,
      globalStress.xy,
      globalStress.yx,
      globalStress.yy,
      matStress.xx,
      matStress.xy,
      matStress.yx,
      matStress.yy,
      sacStress.xx,
      sacStress.xy,
      sacStress.yx,
      sacStress.yy};
}

auto networkV4::protocols::quasiStaticStrainDouble::getCounts(
    const network& _network)
    -> std::tuple<std::size_t, std::size_t, std::size_t>
{
  std::size_t bondCount = 0;
  std::size_t sacrificialCount = 0;
  std::size_t matrixCount = 0;

  auto sacTag = _network.getTags().get("sacrificial");

  for (const auto& [bond, type] : ranges::views::zip(
           _network.getBonds().getBonds(), _network.getBonds().getTypes()))
  {
    if (std::holds_alternative<Forces::VirtualBond>(type)) {
      continue;
    }
    bondCount++;
    if (_network.getBondTags().hasTag(bond.index, sacTag)) {
      sacrificialCount++;
    } else {
      matrixCount++;
    }
  }
  return {bondCount, sacrificialCount, matrixCount};
}

auto networkV4::protocols::quasiStaticStrainDouble::breakData(
    const network& _network) -> std::tuple<double, size_t>
{
  double maxThres = -1e10;
  size_t broken = 0;

  const auto& nodes = _network.getNodes();
  const auto& bonds = _network.getBonds();
  const auto& box = _network.getBox();

  for (const auto& [bond, type, brk] : ranges::views::zip(
           bonds.getBonds(), bonds.getTypes(), bonds.getBreaks()))
  {
    const auto& pos1 = nodes.positions()[bond.src];
    const auto& pos2 = nodes.positions()[bond.dst];
    const auto dist = box.minDist(pos1, pos2);

    if (bonded::visitBreak(brk, dist)) {
      broken++;
    }
    maxThres =
        std::max(maxThres, bonded::visitThreshold(brk, dist).value_or(-1e10));
  }
  return {maxThres, broken};
}

auto networkV4::protocols::quasiStaticStrainDouble::xdoty(
    const std::vector<Utils::vec2d>& _x,
    const std::vector<Utils::vec2d>& _y) -> double
{
  return std::inner_product(_x.begin(),
                            _x.end(),
                            _y.begin(),
                            0.0,
                            std::plus<double>(),
                            [](Utils::vec2d _a, Utils::vec2d _b)
                            { return _a.dot(_b); });
}

//-------------------------------------------------------------------------------------------------

auto networkV4::protocols::quasiStaticStrainDoubleReader::read(
    const toml::value& _config,
    const network& _network,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
    std::shared_ptr<IO::networkDumps::networkDump>& _networkOut)
    -> std::shared_ptr<protocolBase>
{
  auto quasiConfig = toml::find(_config, "QuasiStaticStrain");

  if (!quasiConfig.contains("MaxStrain")) {
    throw std::runtime_error("MaxStrain not specified");
  }
  const double maxStrain = toml::find<double>(quasiConfig, "MaxStrain");

  std::shared_ptr<deform::deformBase> deform = readDeform(quasiConfig);

  integration::AdaptiveParams adaptiveParams = readAdaptive(_config);
  minimisation::minimiserParams minimiserParams = readMinimiser(_config);

  const double tol =
      toml::find_or<double>(quasiConfig, "Tol", config::rootMethods::targetTol);

  const bool errorOnNotSingleBreak = toml::find_or<bool>(
      quasiConfig,
      "ErrorOnNotSingleBreak",
      config::protocols::quasiStaticStrain::errorOnNotSingleBreak);
  const double maxStep =
      toml::find_or<double>(quasiConfig, "MaxStep", maxStrain);
  const bool saveBreaks = toml::find_or<bool>(quasiConfig, "SaveBreaks", false);

  return std::make_shared<quasiStaticStrainDouble>(deform,
                                                   _dataOut,
                                                   _bondsOut,
                                                   _networkOut,
                                                   _network,
                                                   maxStrain,
                                                   tol,
                                                   adaptiveParams,
                                                   minimiserParams,
                                                   errorOnNotSingleBreak,
                                                   maxStep,
                                                   saveBreaks);
}