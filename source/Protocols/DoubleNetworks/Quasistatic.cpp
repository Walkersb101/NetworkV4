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
    networkSavePoints _savePoints)
    : protocolBase(_deform, _dataOut, _bondsOut, _networkOut, _network)
    , m_maxStrain(_maxStrain)
    , m_rootTol(_rootTol)
    , m_params(_params)
    , m_minParams(_minParams)
    , m_errorOnNotSingleBreak(_errorOnNotSingleBreak)
    , m_maxStep(_maxStep)
    , m_savePoints(_savePoints)
    , m_breakMinimiser(*this)
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
      "time",
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
      "ShearStrain",
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
  _network = evalStrain(_network, 0.0);
  logData(_network, "Initial", 0, 0.0, true);

  while (true) {
    auto state = findNextBreak(_network);
    switch (state) {
      case nextBreakState::FoundSingleBreak:
        break;
      case nextBreakState::FoundMultipleBreaks: {
        if (m_errorOnNotSingleBreak)
          throw std::runtime_error(
              "More than one bond above threshold and "
              "errorOnNotSingleBreak is true");
        else
          std::cout << "More than one bond above threshold" << std::endl;
        break;
      }
      case nextBreakState::MaxStrainReached: {
        std::cout << "Max Strain Reached" << std::endl;
        return;
      }
    }

    if (m_savePoints.time) {
      std::visit([&](auto& info) { info.reset(); }, m_savePoints.time.value());
    };

    m_breakMinimiser.minimise(_network);
    if (m_oneBreak)
      break;
  }
}

auto networkV4::protocols::quasiStaticStrainDouble::evalStrain(
    const network& _network, double _targetStrain) -> network
{
  network result = _network;

  const double step = _targetStrain - m_deform->getStrain(result);
  m_deform->strain(result, step);

  minimisation::fire2 minimizer(m_minParams);
  // minimisation::SD minimizer(m_minParams);
  minimizer.minimise(result);
  result.computeForces<false, true>();
  return result;
}

auto networkV4::protocols::quasiStaticStrainDouble::converge(network& _network,
                                                             double _a,
                                                             double _b,
                                                             double _fa,
                                                             double _fb,
                                                             double _tol)
    -> roots::rootState
{
  auto networkA = _network;
  auto networkB = _network;

  if (_a > _b)
    return roots::rootState::MinLargerThanMax;

  if (_fa * _fb > 0.0)
    return roots::rootState::RootNotBracketed;

  roots::ITP solver(_a, _b, _tol);

  for (std::size_t iters = 0; iters < solver.nMax(); ++iters) {
    auto xITP = solver.guessRoot(_a, _b, _fa, _fb);
    if (xITP <= _a || xITP >= _b) {
      _network = networkA;
      std::cout << "Root not bracketed: accepting lower bound" << std::endl;
      return roots::rootState::GuessNotInBracket;
    }

    auto networkITP = evalStrain(_network, xITP);
    auto [maxDistAboveITP, breakCountITP] = breakData(networkITP);
    if (breakCountITP > 0) {
      _b = xITP;
      _fb = maxDistAboveITP;
      networkB = networkITP;
    } else {
      _a = xITP;
      _fa = maxDistAboveITP;
      networkA = networkITP;
    }

    if (std::abs(_b - _a) < 2 * _tol) {
      _network = networkB;
      return roots::rootState::converged;
    }
  }
  std::cout << "Max iterations reached: ";
  if (_fb > 0.0) {
    std::cout << "fb > 0 accepting upper bound" << std::endl;
    _network = networkB;
  } else {
    std::cout << "fa > 0 accepting lower bound" << std::endl;
    _network = networkA;
  }
  return roots::rootState::MaxIterationsReached;
}

auto networkV4::protocols::quasiStaticStrainDouble::findNextBreak(
    network& _network) -> nextBreakState
{
  bool breakFound = false;
  while (true) {
    double a = m_deform->getStrain(_network);
    double b = std::min(a + m_maxStep, m_maxStrain);

    if (a >= m_maxStrain - 1e-10)
      return nextBreakState::MaxStrainReached;

    auto [fa, breakCountA] = breakData(_network);
    if (breakCountA == 1)
      return nextBreakState::FoundSingleBreak;
    if (breakCountA > 1)
      return nextBreakState::FoundMultipleBreaks;

    m_strainCount++;
    auto bNetwork = evalStrain(_network, b);
    auto [fb, breakCountB] = breakData(bNetwork);

    if (breakCountB == 0) {
      _network = bNetwork;
      logData(_network, "Strain", 0, 0.0, true);
      continue;
    }

    auto state = converge(_network, a, b, fa, fb, m_rootTol);
    switch (state) {
      case roots::rootState::MinLargerThanMax:
        throw std::runtime_error("Min larger than max");
      case roots::rootState::RootNotBracketed:
        throw std::runtime_error("Root not bracketed");
    }
    auto [maxDistAbove, breakCount] = breakData(_network);
    if (breakCount == 1)
      return nextBreakState::FoundSingleBreak;
    if (breakCount > 1)
      return nextBreakState::FoundMultipleBreaks;

    std::cout << "Converged with zero breaks retrying with lower bound"
              << std::endl;
    logData(_network, "LowerBound", 0, 0.0, false);
  }
}

auto networkV4::protocols::quasiStaticStrainDouble::genTimeData(
    const network& _network,
    const std::string& _reason,
    size_t _breakCount,
    double _t) -> std::vector<IO::timeSeries::writeableTypes>
{
  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress =
      stresses.getFirst(_network.getTags().get("sacrificial"));
  const auto& matStress = stresses.getFirst(_network.getTags().get("matrix"));

  const auto& box = _network.getBox();
  auto counts = getCounts(_network);

  return {_reason,
          m_strainCount,
          _breakCount,
          _t,
          box.getLx(),
          box.getLy(),
          box.shearStrain(),
          _network.getElongationStrain()[0],
          _network.getElongationStrain()[1],
          std::get<0>(counts),
          std::get<1>(counts),
          std::get<2>(counts),
          globalStress(0, 0),
          globalStress(0, 1),
          globalStress(1, 0),
          globalStress(1, 1),
          matStress(0, 0),
          matStress(0, 1),
          matStress(1, 0),
          matStress(1, 1),
          sacStress(0, 0),
          sacStress(0, 1),
          sacStress(1, 0),
          sacStress(1, 1)};
}

auto networkV4::protocols::quasiStaticStrainDouble::genBondData(
    const network& _network, const breakInfo& _bond, double _t)
    -> std::vector<IO::timeSeries::writeableTypes>
{
  auto sacTag = _network.getTags().get("sacrificial");
  auto matTag = _network.getTags().get("matrix");

  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress = stresses.getFirst(sacTag);
  const auto& matStress = stresses.getFirst(matTag);

  const auto& box = _network.getBox();
  auto counts = getCounts(_network);

  const auto& nodes = _network.getNodes();
  const auto bonds = _network.getBonds();

  const auto& binfo = std::get<0>(_bond);
  const auto& type = std::get<1>(_bond);
  const auto& brk = std::get<2>(_bond);
  const auto tags = std::get<3>(_bond);

  bool harmonic = std::holds_alternative<Forces::HarmonicBond>(type);
  bool strainBreak = std::holds_alternative<BreakTypes::StrainBreak>(brk);
  bool sacrificial = Utils::Tags::hasTag(tags, sacTag);

  const auto& pos1 = nodes.positions()[binfo.src];
  const auto& pos2 = nodes.positions()[binfo.dst];
  size_t bondSrc = nodes.indices()[binfo.src];
  size_t bondDst = nodes.indices()[binfo.dst];

  double r = box.minDist(pos1, pos2).norm();
  double r0 = harmonic ? std::get<Forces::HarmonicBond>(type).r0() : 0.0;

  return {
      m_strainCount,
      _t,
      sacrificial ? "Sacrificial" : "Matrix",
      harmonic ? std::get<Forces::HarmonicBond>(std::get<1>(_bond)).k() : 0.0,
      strainBreak ? std::get<BreakTypes::StrainBreak>(brk).lambda() : 0.0,
      harmonic ? std::get<Forces::HarmonicBond>(std::get<1>(_bond)).r0() : 0.0,
      harmonic ? (r - r0) / r0 : 0.0,
      bondSrc,
      bondDst,
      pos1[0],
      pos1[1],
      pos2[0],
      pos2[1],
      box.getLx(),
      box.getLy(),
      box.shearStrain(),
      _network.getElongationStrain()[0],
      _network.getElongationStrain()[1],
      std::get<0>(counts),
      std::get<1>(counts),
      std::get<2>(counts),
      globalStress(0, 0),
      globalStress(0, 1),
      globalStress(1, 0),
      globalStress(1, 1),
      matStress(0, 0),
      matStress(0, 1),
      matStress(1, 0),
      matStress(1, 1),
      sacStress(0, 0),
      sacStress(0, 1),
      sacStress(1, 0),
      sacStress(1, 1)};
}

void networkV4::protocols::quasiStaticStrainDouble::logData(
    network& _network,
    const std::string& _reason,
    double _breakcount,
    double _t,
    bool _writeDump)
{
  _network.computeForces<false, true>();
  m_dataOut->write(genTimeData(_network, _reason, _breakcount, _t));
  if (_writeDump) {
    m_networkOut->save(_network, m_strainCount, _t, _reason);
  }
}

auto networkV4::protocols::quasiStaticStrainDouble::getCounts(
    const network& _network)
    -> std::tuple<std::size_t, std::size_t, std::size_t>
{
  std::size_t bondCount = 0;
  std::size_t sacrificialCount = 0;
  std::size_t matrixCount = 0;

  auto sacTag = _network.getTags().get("sacrificial");
  const auto& bonds = _network.getBonds();

  for (const auto& [bond, type, tags] :
       ranges::views::zip(bonds.getBonds(), bonds.getTypes(), bonds.getTags()))
  {
    if (std::holds_alternative<Forces::VirtualBond>(type)) {
      continue;
    }
    bondCount++;
    if (Utils::Tags::hasTag(tags, sacTag)) {
      sacrificialCount++;
    } else {
      matrixCount++;
    }
  }
  return {bondCount, sacrificialCount, matrixCount};
}

auto networkV4::protocols::quasiStaticStrainDouble::checkIfNeedToSave(
    const network& _network, double _t) -> std::optional<std::string>
{
  std::optional<std::string> reason = std::nullopt;
  if (m_deform->getStrain(_network) >= m_savePoints.strain)
    reason = "Strained";
  if (_t >= m_savePoints.time)
    reason = "Time";
  if (Utils::brokenBonds(_network) >= m_savePoints.breakCount)
    reason = "Broken";
  if (m_strainCount >= m_savePoints.strainStep)
    reason = "Strain";

  return reason;
}

auto networkV4::protocols::quasiStaticStrainDouble::processBreakQueue(
    network& _network, double _t) -> size_t
{
  if (_network.getBreakQueue().empty())
    return 0;
  _network.computeForces<false, true>();
  size_t breakCount = _network.getBreakQueue().size();
  while (!_network.getBreakQueue().empty()) {
    const auto& broken = _network.getBreakQueue().front();
    m_bondsOut->write(genBondData(_network, broken, _t));
    _network.getBreakQueue().pop_front();
  }
  return breakCount;
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

//-------------------------------------------------------------------------------------------------

void networkV4::protocols::quasiStaticStrainDouble::relaxBreak::minimise(
    network& _network)
{
  _network.computeForces<true, true>();
  size_t breakCount = m_protocol.processBreakQueue(_network, 0.0);

  bool dump = m_protocol.checkIfNeedToSave(_network).has_value();
  m_protocol.logData(_network, "Start", breakCount, 0.0, dump);

  double Ecurr = _network.getEnergy();
  double Eprev = Ecurr;

  integration::AdaptiveOverdampedEulerHeun stepper(1.0, m_params);

  double t = 0.0;
  size_t iter = 0;
  while (iter++ < m_maxIter) {
    Eprev = Ecurr;

    auto status = hybridStep(_network, stepper);
    if (!status) {
      switch (status.error()) {
        case lineSearch::lineSearchState::DirectionNotDescent:
          throw std::runtime_error("Direction not descent");
        case lineSearch::lineSearchState::zeroquad:
          throw std::runtime_error("Zero quad");
        case lineSearch::lineSearchState::zeroAlpha:
          throw std::runtime_error("Zero alpha");
        case lineSearch::lineSearchState::zeroforce:
          std::cout << "Zero force" << std::endl;
          break;
        default:
          throw std::runtime_error("Error in Line Search");
      }
      break;
    }

    _network.computeForces<true, false>();
    Ecurr = _network.getEnergy();
    t += status.value();

    bool brokenInStep = !_network.getBreakQueue().empty();
    breakCount += m_protocol.processBreakQueue(_network, t);

    auto reason = m_protocol.checkIfNeedToSave(_network);
    if (reason) {
      m_protocol.logData(_network, reason.value(), breakCount, t, true);
    }
    double fdotf = Utils::Math::xdoty(_network.getNodes().forces(),
                                      _network.getNodes().forces());
    if (!brokenInStep && converged(fdotf, Ecurr, Eprev)) {
      break;
    }
  }
  m_protocol.logData(_network, "End", breakCount, t, dump);
}

auto networkV4::protocols::quasiStaticStrainDouble::relaxBreak::hybridStep(
    network& _network, auto& _stepper)
    -> tl::expected<double, lineSearch::lineSearchState>
{
  auto rk = _network.getNodes().positions();
  auto fk = _network.getNodes().forces();
  double Eoriginal = _network.getEnergy();

  _stepper.step(_network);
  _network.computeForces();
  double Ecurr = _network.getEnergy();
  if (Ecurr < Eoriginal)
    return _stepper.getDt();

  _network.getNodes().positions() = rk;
  _network.getNodes().forces() = fk;

  const double zeta = _stepper.getZeta();  // TODO: allow changing zeta
  lineSearch::lineSearchQuad lineSearch(zeta * m_params.dtMax);
  auto status = lineSearch.search(fk, _network);
  _network.computeForces();
  if (status)
    return status.value() / zeta;
  else
    return tl::make_unexpected(status.error());
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

  auto saveConfig = readSavePoints(quasiConfig);

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
                                                   saveConfig);
}

auto networkV4::protocols::quasiStaticStrainDoubleReader::readSavePoints(
    const toml::value& _config) -> quasiStaticStrainDouble::networkSavePoints
{
  auto checkConfig =
      [](const toml::value& saveConfig, const std::string& baseName)
  {
    return saveConfig.contains(baseName + "Scale")
        || saveConfig.contains(baseName + "Step")
        || saveConfig.contains(baseName + "Start");
  };

  auto getScale = [](const toml::value& saveConfig, const std::string& baseName)
  {
    return toml::find_or<std::string>(saveConfig, baseName + "Scale", "Linear")
            == "Linear"
        ? Utils::Output::stepScale::Linear
        : Utils::Output::stepScale::Log;
  };

  quasiStaticStrainDouble::networkSavePoints savePoints;
  if (_config.contains("SaveNetworkOn")) {
    const auto& saveConfig = toml::find(_config, "SaveNetworkOn");

    if (toml::find_or<bool>(saveConfig, "StrainCount", false)) {
      savePoints.strainStep = Utils::Output::linearStepInfo<size_t>(0, 1);
    }
    if (checkConfig(saveConfig, "StrainCount")) {
      auto scale = getScale(saveConfig, "StrainCount");
      auto start = toml::find_or<size_t>(saveConfig, "StrainCountStart", 0);
      switch (scale) {
        case Utils::Output::stepScale::Linear: {
          auto step = toml::find_or<size_t>(saveConfig, "StrainCountStep", 1);
          savePoints.strainStep =
              Utils::Output::linearStepInfo<size_t>(start, step);
          break;
        }
        case Utils::Output::stepScale::Log: {
          auto step = toml::find_or<double>(saveConfig, "StrainCountStep", 2.0);
          savePoints.strainStep =
              Utils::Output::logStepInfo<size_t, double>(start, step);
          break;
        }
        default:
          break;
      }
    }

    if (toml::find_or<bool>(saveConfig, "BreakCount", false)) {
      savePoints.breakCount = Utils::Output::linearStepInfo<size_t>(1, 1);
    }
    if (checkConfig(saveConfig, "BreakCount")) {
      auto scale = getScale(saveConfig, "BreakCount");
      auto start = toml::find_or<size_t>(saveConfig, "BreakCountStart", 1);
      switch (scale) {
        case Utils::Output::stepScale::Linear: {
          auto step = toml::find_or<size_t>(saveConfig, "BreakCountStep", 1);
          savePoints.breakCount =
              Utils::Output::linearStepInfo<size_t>(start, step);
          break;
        }
        case Utils::Output::stepScale::Log: {
          auto step = toml::find_or<double>(saveConfig, "BreakCountStep", 2.0);
          savePoints.breakCount =
              Utils::Output::logStepInfo<size_t, double>(start, step);
          break;
        }
        default:
          break;
      }
    }

    if (toml::find_or<bool>(saveConfig, "Time", false)) {
      savePoints.time = Utils::Output::linearStepInfo<double>(0.0, 2.0);
    }
    if (checkConfig(saveConfig, "Time")) {
      auto scale = getScale(saveConfig, "Time");
      auto start = toml::find_or<double>(saveConfig, "TimeStart", 0.0);
      switch (scale) {
        case Utils::Output::stepScale::Linear: {
          auto step = toml::find_or<double>(saveConfig, "TimeStep", 2.0);
          savePoints.time = Utils::Output::linearStepInfo<double>(start, step);
          break;
        }
        case Utils::Output::stepScale::Log: {
          auto step = toml::find_or<double>(saveConfig, "TimeStep", 1.1);
          savePoints.time =
              Utils::Output::logStepInfo<double, double>(start, step);
          break;
        }
        default:
          break;
      }
    }

    if (toml::find_or<bool>(saveConfig, "Strain", false)) {
      savePoints.strain = Utils::Output::linearStepInfo<double>(0.0, 0.1);
    }
    if (checkConfig(saveConfig, "Strain")) {
      auto scale = getScale(saveConfig, "Strain");
      auto start = toml::find_or<double>(saveConfig, "StrainStart", 0.0);
      switch (scale) {
        case Utils::Output::stepScale::Linear: {
          auto step = toml::find_or<double>(saveConfig, "StrainStep", 0.1);
          savePoints.strain =
              Utils::Output::linearStepInfo<double>(start, step);
          break;
        }
        case Utils::Output::stepScale::Log: {
          auto step = toml::find_or<double>(saveConfig, "StrainStep", 1.1);
          savePoints.strain =
              Utils::Output::logStepInfo<double, double>(start, step);
          break;
        }
        default:
          break;
      }
    }
  }
  return savePoints;
}