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
  m_dataOut->write(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  bool stop = false;

  while (!stop) {
    auto newNetwork = findNextBreak(_network, m_errorOnNotSingleBreak);
    if (!newNetwork) {
      switch (newNetwork.error()) {
        case nextBreakErrors::BreakAtLowerBound:
          throw std::runtime_error("Break at lower bound");
        case nextBreakErrors::DidNotConverge:
          throw std::runtime_error("Did not converge");
        case nextBreakErrors::convergedWithMoreThanOneBreak:
          throw std::runtime_error(
              "More than one bond above threshold and "
              "errorOnNotSingleBreak is true");
        case nextBreakErrors::ConvergedWithZeroBreaks: {
          // TODO: Log
        }
        case nextBreakErrors::MaxStrainReached: {
          stop = true;
          continue;
        }
        default:
          throw std::runtime_error("Unknown error");
      }
    }

    m_t = 0.0;
    if (m_savePoints.time) {
      std::visit([&](auto& info) { info.reset(); }, m_savePoints.time.value());
    }

    m_strainCount++;
    _network = newNetwork.value();
    _network.computeForces<true, true>();
    auto breakCount = processBreakQueue(_network);
    auto reason = checkIfNeedToSave(_network);
    m_dataOut->write(genTimeData(_network, "Start", breakCount));
    if (reason) {
      m_networkOut->save(_network,
                         m_strainCount,
                         0.0,
                         "Start-" + std::to_string(m_strainCount));
    }

    breakCount = relaxBreak(_network, breakCount);
    m_dataOut->write(genTimeData(_network, "End", breakCount));
    if (reason) {
      m_networkOut->save(
          _network, m_strainCount, m_t, "End-" + std::to_string(m_strainCount));
    }

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

auto networkV4::protocols::quasiStaticStrainDouble::converge(
    const network& _aNetwork,
    const network& _bNetwork,
    double _a,
    double _b,
    double _fa,
    double _fb,
    double _tol) -> tl::expected<network, roots::rootErrors>
{
  auto networkB = _bNetwork;

  if (_a > _b)
    return tl::make_unexpected(roots::rootErrors::MinLargerThanMax);

  if (_fa * _fb > 0.0)
    return tl::make_unexpected(roots::rootErrors::RootNotBracketed);

  roots::ITP solver(_a, _b, _tol);

  for (std::size_t iters = 0; iters < solver.nMax(); ++iters) {
    auto xITP = solver.guessRoot(_a, _b, _fa, _fb);

    auto networkITP = evalStrain(_aNetwork, xITP);
    auto [maxDistAboveITP, breakCountITP] = breakData(networkITP);
    if (maxDistAboveITP >= 0.) {
      _b = xITP;
      _fb = maxDistAboveITP;
      networkB = std::move(networkITP);
    } else {
      _a = xITP;
      _fa = maxDistAboveITP;
    }

    if (std::abs(_b - _a) < 2 * _tol) {
      return networkB;
    }
  }
  return tl::make_unexpected(roots::rootErrors::MaxIterationsReached);
}

auto networkV4::protocols::quasiStaticStrainDouble::findNextBreak(
    const network& _network, bool _singleBreak)
    -> tl::expected<network, nextBreakErrors>
{
  auto startStrain = m_deform->getStrain(_network);
  auto resultNetwork = evalStrain(_network, startStrain);

  if (get<size_t>(breakData(resultNetwork)) > 0)
    return tl::make_unexpected(nextBreakErrors::BreakAtLowerBound);

  bool breakFound = false;
  while (!breakFound) {
    double a = m_deform->getStrain(resultNetwork);
    double b = std::min(a + m_maxStep, m_maxStrain);

    if (a == m_maxStrain)
      return tl::make_unexpected(nextBreakErrors::MaxStrainReached);

    auto [fa, breakCountA] = breakData(resultNetwork);

    if (breakCountA > 0)
      return tl::make_unexpected(nextBreakErrors::BreakAtLowerBound);

    auto bNetwork = evalStrain(resultNetwork, b);
    auto [fb, breakCountB] = breakData(bNetwork);

    if (breakCountB == 0) {
      resultNetwork = std::move(bNetwork);
      m_strainCount++;
      m_dataOut->write(genTimeData(resultNetwork, "Strain", 0));
      auto reason = checkIfNeedToSave(resultNetwork);
      if (reason) {
        m_networkOut->save(resultNetwork, m_strainCount, 0.0, reason.value());
      }
      continue;
    }

    auto converged = converge(resultNetwork, bNetwork, a, b, fa, fb, m_rootTol);
    if (!converged) {
      return tl::make_unexpected(nextBreakErrors::DidNotConverge);
    }

    resultNetwork = converged.value();
    auto [maxDistAbove, breakCount] = breakData(resultNetwork);

    if (breakCount == 0) {
      return tl::make_unexpected(nextBreakErrors::ConvergedWithZeroBreaks);
    }
    if (_singleBreak && breakCount > 1) {
      return tl::make_unexpected(
          nextBreakErrors::convergedWithMoreThanOneBreak);
    }

    breakFound = true;
  }
  return resultNetwork;
}

auto networkV4::protocols::quasiStaticStrainDouble::relaxBreak(
    network& _network, size_t startBreaks) -> size_t
{
  double Ecurr = _network.getEnergy();
  double Eprev = Ecurr;
  size_t breakCount = startBreaks;

  integration::AdaptiveOverdampedEulerHeun stepper(1.0, m_params);

  size_t iter = 0;
  while (iter++ < m_minParams.maxIter) {
    Eprev = Ecurr;

    auto status = hybridStep(_network, stepper);
    double dt = status.value_or(0.0);
    if (!status) {
      if (status.error() == lineSearch::lineSearchState::zeroAlpha
          || status.error() == lineSearch::lineSearchState::zeroquad)
      {
        break;
      } else {
        throw std::runtime_error("lineSearch error");
      }
    }

    _network.computeForces<true, false>();
    Ecurr = _network.getEnergy();
    m_t += dt;

    bool brokenInStep = !_network.getBreakQueue().empty();
    breakCount += processBreakQueue(_network);
    auto reason = checkIfNeedToSave(_network);
    if (reason) {
      m_dataOut->write(genTimeData(_network, reason.value(), breakCount));
      m_networkOut->save(_network, m_strainCount, m_t, reason.value());
    }

    if (!brokenInStep
        && fabs(Ecurr - Eprev) < m_minParams.Etol * 0.5
                * (fabs(Ecurr) + fabs(Eprev) + minimisation::EPS_ENERGY))
    {
      break;
    }

    double fdotf = Utils::Math::xdoty(_network.getNodes().forces(),
                                      _network.getNodes().forces());
    if (!brokenInStep && fdotf < m_minParams.Ftol * m_minParams.Ftol) {
      break;
    }
    // std::cout << std::setprecision(10) << iter << " " << error << " "
    //           << _network.getEnergy() << " " << integrator.getDt() << " "
    //           << tools::maxLength(_network.getNodes().forces()) << " "
    //           << tools::maxAbsComponent(_network.getNodes().forces())
    //           << std::endl;
  }
  _network.computeForces<false, true>();
  return breakCount;
}

auto networkV4::protocols::quasiStaticStrainDouble::genTimeData(
    const network& _network, const std::string& _reason, size_t _breakCount)
    -> std::vector<IO::timeSeries::writeableTypes>
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
          m_t,
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
    const network& _network, const breakInfo& _bond)
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
      m_t,
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
    const network& _network) -> std::optional<std::string>
{
  std::optional<std::string> reason = std::nullopt;
  if (m_deform->getStrain(_network) >= m_savePoints.strain)
    reason = "Strained";
  if (m_t >= m_savePoints.time)
    reason = "Time";
  if (Utils::brokenBonds(_network) >= m_savePoints.breakCount)
    reason = "Broken";
  if (m_strainCount >= m_savePoints.strainStep)
    reason = "Strain";

  return reason;
}

auto networkV4::protocols::quasiStaticStrainDouble::processBreakQueue(
    network& _network) -> size_t
{
  size_t breakCount = _network.getBreakQueue().size();
  while (!_network.getBreakQueue().empty()) {
    const auto& broken = _network.getBreakQueue().front();
    m_bondsOut->write(genBondData(_network, broken));
    _network.getBreakQueue().pop_front();
  }
  return breakCount;
}

auto networkV4::protocols::quasiStaticStrainDouble::hybridStep(
    network& _network, auto& _stepper)
    -> tl::expected<double, lineSearch::lineSearchState>
{
  auto rk = _network.getNodes().positions();
  double Eoriginal = _network.getEnergy();

  _stepper.step(_network);
  _network.computeForces();
  double Ecurr = _network.getEnergy();
  if (Ecurr < Eoriginal)
    return _stepper.getDt();

  _network.getNodes().positions() = rk;
  _network.computeForces();

  const double zeta = _stepper.getZeta();  // TODO: allow changing zeta
  lineSearch::lineSearchQuad lineSearch(zeta * m_params.dtMax);
  auto h = _network.getNodes().forces();
  auto status = lineSearch.search(h, _network);
  _network.computeForces();
  return status;
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