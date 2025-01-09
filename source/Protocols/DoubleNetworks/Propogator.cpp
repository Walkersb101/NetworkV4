#include "Propogator.hpp"
#include <iostream>

networkV4::protocols::propogatorDouble::propogatorDouble(
    std::shared_ptr<deform::deformBase>& _deform,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
    std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
    const network& _network,
    const std::vector<double>& _strains,
    double _rootTol,
    integration::AdaptiveParams _params,
    const minimisation::minimiserParams& _minParams,
    double _maxStep)
    : protocolBase(_deform, _dataOut, _bondsOut, _networkOut, _network)
    , m_strains(_strains)
    , m_rootTol(_rootTol)
    , m_params(_params)
    , m_minParams(_minParams)
    , m_maxStep(_maxStep)
{
  std::vector<IO::timeSeries::writeableTypes> dataHeader = {
      "Reason",
      "StrainCount",
      "Domainx",
      "Domainy",
      "ShearStrain",
      "ElongationStrainX",
      "ElongationStrainY",
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
      "Type",
      "mu",
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
      "ElongationStrainY"};

  _dataOut->write(dataHeader);
  _bondsOut->write(bondHeader);
}

void networkV4::protocols::propogatorDouble::run(network& _network)
{
  relax(_network);
  m_dataOut->write(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  if (m_strains.empty()) {
    runLambda(_network);
  } else {
    runStrain(_network);
  }
}

void networkV4::protocols::propogatorDouble::runLambda(network& _network)
{
  bool foundBreak = false;
  while (!foundBreak) {
    foundBreak = findSingleBreak(_network);
    std::cout << m_deform->getStrain(_network) << std::endl;
  }
  m_strainCount = 1;

  _network.computeBreaks();
  while (!_network.getBreakQueue().empty()) {
    const auto& broken = _network.getBreakQueue().front();
    m_bondsOut->write(genBondData(_network, broken));
    _network.getBreakQueue().pop_front();
  }

  m_dataOut->write(genTimeData(_network, "Start", 1));
  m_networkOut->save(_network, 1, 0.0, "Start");
  relax(_network);
  m_dataOut->write(genTimeData(_network, "End", 1));
  m_networkOut->save(_network, 1, 1.0, "End");
}

void networkV4::protocols::propogatorDouble::runStrain(network& _network)
{
  std::sort(m_strains.begin(), m_strains.end());
  network SavedNetwork = _network;

  for (const auto& targetStrain : m_strains) {
    while (m_deform->getStrain(_network) <= targetStrain - 1e-14) {
      const double subStepStrain =
          std::min(targetStrain, m_deform->getStrain(_network) + m_maxStep);
      evalStrain(_network, subStepStrain);
      std::cout << m_deform->getStrain(_network) << std::endl;
    }
    m_strainCount++;
    SavedNetwork = _network;

    breakMostStrained(_network, _network.getTags().get("sacrificial"));

    m_dataOut->write(genTimeData(_network, "Start", 1));
    m_networkOut->save(_network, m_strainCount, 0.0, "Start");

    relax(_network);

    m_dataOut->write(genTimeData(_network, "End", 1));
    m_networkOut->save(_network, m_strainCount, 1.0, "End");

    _network = SavedNetwork;
  }
}

void networkV4::protocols::propogatorDouble::evalStrain(network& _network,
                                                        double _targetStrain)
{
  const double step = _targetStrain - m_deform->getStrain(_network);
  m_deform->strain(_network, step);
  relax(_network);
}

void networkV4::protocols::propogatorDouble::relax(network& _network)
{
  //minimisation::AdaptiveHeunDecent minimizer(m_minParams, m_params);
  minimisation::fire2 minimizer(m_minParams);
  minimizer.minimise(_network);
}

auto networkV4::protocols::propogatorDouble::getMaxDataIndex(
    network& _network, const Utils::Tags::tagFlags& _filter) -> size_t
{
  auto filter = [&_filter](const Utils::Tags::tagFlags& _tags) -> bool
  { return Utils::Tags::hasTagAny(_tags, _filter); };

  const auto& bonds = _network.getBonds();
  const auto& positions = _network.getNodes().positions();
  const auto& box = _network.getBox();

  double max = -1e10;
  size_t maxIndex = 0;
  for (const auto [i, bond, type, brk, tags] :
       ranges::views::zip(ranges::views::iota(0ul, bonds.size()),
                          bonds.getBonds(),
                          bonds.getTypes(),
                          bonds.getBreaks(),
                          bonds.getTags()))
  {
    if (filter(tags)) {
      const auto& pos1 = positions[bond.src];
      const auto& pos2 = positions[bond.dst];
      const auto dist = box.minDist(pos1, pos2);

      auto val = bonded::visitData(brk, dist);
      if (val && val.value() > max) {
        max = val.value();
        maxIndex = i;
      }
    }
  }
  return maxIndex;
}

void networkV4::protocols::propogatorDouble::breakMostStrained(
    network& _network, const Utils::Tags::tagFlags& _filter)
{
  size_t maxIndex = getMaxDataIndex(_network, _filter);
  auto& bonds = _network.getBonds();
  auto& binfo = bonds.getBonds()[maxIndex];
  auto& type = bonds.getTypes()[maxIndex];
  auto& brk = bonds.getBreaks()[maxIndex];
  auto& tags = bonds.getTags()[maxIndex];

  breakInfo b(binfo, type, brk, tags);
  m_bondsOut->write(genBondData(_network, b));

  type = Forces::VirtualBond {};
  brk = BreakTypes::None {};
  tags.set(BROKEN_TAG_INDEX);
}

auto networkV4::protocols::propogatorDouble::breakData(const network& _network)
    -> std::tuple<double, size_t>
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

auto networkV4::protocols::propogatorDouble::findSingleBreak(network& _network)
    -> bool
{
  network testNetwork = _network;
  double maxDistAboveA, maxDistAboveB;
  size_t breakCountA, breakCountB;
  bool converged = false;

  double a = m_deform->getStrain(_network);
  double b = a + m_maxStep;

  testNetwork = _network;
  evalStrain(testNetwork, a);
  std::tie(maxDistAboveA, breakCountA) = breakData(testNetwork);

  // TODO : add Log
  if (maxDistAboveA > 0.0) {
    return false;
  }

  testNetwork = _network;
  evalStrain(testNetwork, b);
  std::tie(maxDistAboveB, breakCountB) = breakData(testNetwork);

  // TODO : add Log
  // Step not large enough to break
  if (maxDistAboveB < 0.0) {
    _network = testNetwork;
    return false;
  }

  if (maxDistAboveA * maxDistAboveB > 0.0) {
    throw std::invalid_argument("Root not bracketed");
    return false;
  }

  roots::ITP solver(a, b, m_rootTol);
  for (size_t iters = 0; iters < solver.nMax(); ++iters) {
    double xITP = solver.guessRoot(a, b, maxDistAboveA, maxDistAboveB);

    testNetwork = _network;
    evalStrain(testNetwork, xITP);
    auto [maxDistAboveITP, breakCountITP] = breakData(testNetwork);
    if (maxDistAboveITP >= 0.) {
      b = xITP;
      maxDistAboveB = maxDistAboveITP;
    } else {
      a = xITP;
      maxDistAboveA = maxDistAboveITP;
    }
    if (std::abs(b - a) < 2 * m_rootTol) {
      break;
    }
  }
  _network = testNetwork;
  return true;
}

auto networkV4::protocols::propogatorDouble::genTimeData(
    const network& _network,
    const std::string& _reason,
    size_t _breakCount) -> std::vector<IO::timeSeries::writeableTypes>
{
  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress =
      stresses.getFirst(_network.getTags().get("sacrificial"));
  const auto& matStress = stresses.getFirst(_network.getTags().get("matrix"));

  const auto& box = _network.getBox();

  return {_reason,
          m_strainCount,
          _breakCount,
          box.getLx(),
          box.getLy(),
          box.shearStrain(),
          _network.getElongationStrain().x,
          _network.getElongationStrain().y,
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

auto networkV4::protocols::propogatorDouble::genBondData(
    const network& _network,
    const breakInfo& _bond) -> std::vector<IO::timeSeries::writeableTypes>
{
  auto sacTag = _network.getTags().get("sacrificial");
  auto matTag = _network.getTags().get("matrix");

  const auto& stresses = _network.getStresses();
  const auto& globalStress = stresses.total();
  const auto& sacStress = stresses.getFirst(sacTag);
  const auto& matStress = stresses.getFirst(matTag);

  const auto& box = _network.getBox();
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
      box.shearStrain(),
      _network.getElongationStrain().x,
      _network.getElongationStrain().y};
}

// -------------------------------------------------------------------------------------------------

auto networkV4::protocols::propogatorDoubleReader::read(
    const toml::value& _config,
    const network& _network,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
    std::shared_ptr<IO::networkDumps::networkDump>& _networkOut)
    -> std::shared_ptr<protocolBase>
{
  auto propConfig = toml::find(_config, "Propogator");

  std::vector<double> strains =
      toml::find_or<std::vector<double>>(propConfig, "Strains", {});

  std::shared_ptr<deform::deformBase> deform = readDeform(propConfig);

  integration::AdaptiveParams adaptiveParams = readAdaptive(_config);
  minimisation::minimiserParams minimiserParams = readMinimiser(_config);

  const double tol =
      toml::find_or<double>(propConfig, "Tol", config::rootMethods::targetTol);

  const double maxStep =
      toml::find_or<double>(propConfig, "MaxStep", config::protocols::maxStep);

  return std::make_shared<propogatorDouble>(deform,
                                            _dataOut,
                                            _bondsOut,
                                            _networkOut,
                                            _network,
                                            strains,
                                            tol,
                                            adaptiveParams,
                                            minimiserParams,
                                            maxStep);
}