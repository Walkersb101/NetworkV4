#include "Propogator.hpp"

#include "Integration/Minimization.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"

networkV4::propogator::propogator() {}

networkV4::propogator::propogator(std::vector<double> _strains,
                                  StrainType _strainType)
    : m_strains(_strains)
    , protocol(_strainType)
{
}

networkV4::propogator::propogator(const std::vector<double>& _strains,
                                  StrainType _strainType,
                                  bondType _breakType,
                                  double _esp,
                                  double _tol,
                                  double _maxStep,
                                  bool _respectLambda)
    : m_strains(_strains)
    , m_breakType(_breakType)
    , m_esp(_esp)
    , m_tol(_tol)
    , m_maxStep(_maxStep)
    , protocol(_strainType)
    , m_respectLambda(_respectLambda)
{
}

networkV4::propogator::~propogator() {}

void networkV4::propogator::run(network& _network)
{
  relax(_network);
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  if (m_respectLambda) {
    runLambda(_network);
  } else {
    runStrain(_network);
  }
}

void networkV4::propogator::initIO(const network& _network,
                                   std::unique_ptr<dataOut>& _dataOut,
                                   std::unique_ptr<networkOut>& _networkOut)
{
  auto types = tools::uniqueBondTypes(_network.getBonds());

  std::vector<std::string> dataHeader = {
      "Reason",
      "StrainCount",
      "Domainx",
      "Domainy",
      enumString::strainType2Str.at(m_strainType) + "Strain",
  };

  std::vector<std::string> bondHeader = {
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
      enumString::strainType2Str.at(m_strainType) + "Strain",
  };

  tensorData::addBondCountByTypeHeader(dataHeader, types);
  tensorData::addBondCountByTypeHeader(bondHeader, types);

  tensorData::addStressByTypeHeader(bondHeader, types);
  tensorData::addStressByTypeHeader(dataHeader, types);

  _dataOut->initFiles(dataHeader, bondHeader);
  m_dataOut = _dataOut.get();
  m_networkOut = _networkOut.get();
}

void networkV4::propogator::runLambda(network& _network)
{
  bool foundBreak = false;
  network savedNetwork = _network;
  bonds& bonds = _network.getBonds();
  while (!foundBreak) {
    foundBreak = findSingleBreak(_network);
  }
  m_strainCount = 1;
  breakMostStrained(_network, m_breakType);
  m_dataOut->writeTimeData(genTimeData(_network, "Start", 1));
  m_networkOut->save(_network, 1, 0.0, "Start");
  relax(_network);
  m_dataOut->writeTimeData(genTimeData(_network, "End", 1));
  m_networkOut->save(_network, 1, 1.0, "End");
}

void networkV4::propogator::runStrain(network& _network)
{
  std::sort(m_strains.begin(), m_strains.end());
  network SavedNetwork;
  for (const auto& strain : m_strains) {
    m_strainCount++;
    evalStrain(_network, strain);
    SavedNetwork = _network;
    breakMostStrained(_network, m_breakType);
    m_dataOut->writeTimeData(genTimeData(_network, "Start", 1));
    m_networkOut->save(_network, m_strainCount, 0.0, "Start");
    relax(_network);
    m_dataOut->writeTimeData(genTimeData(_network, "End", 1));
    m_networkOut->save(_network, m_strainCount, 1.0, "End");
    _network = SavedNetwork;
  }
}

void networkV4::propogator::evalStrain(network& _network, double _strain)
{
  while (getStrain(_network) < _strain - 1e-15) {
    const double strainStep =
        std::min(_strain - getStrain(_network), m_maxStep);
    strain(_network, strainStep);
    _network.computeForces();
    relax(_network);
  }
}

void networkV4::propogator::relax(network& _network)
{
  FireMinimizer minimizer(m_tol);
  minimizer.integrate(_network);
}

auto networkV4::propogator::getMostStrained(network& _network,
                                            bondType _type) -> size_t
{
  auto isType = [&_type](const bond& b)
  { return (_type == bondType::any) || (b.type() == _type); };
  auto isConnected = [](const bond& b) { return b.connected(); };

  size_t maxIndex = 0;
  double maxStrain = 0.0;
  for (size_t i = 0; i < _network.getBonds().size(); ++i) {
    const auto& b = _network.getBonds().get(i);
    if (!isType(b) || !isConnected(b)) {
      continue;
    }
    const double strain = _network.bondStrain(b);
    if (strain > maxStrain) {
      maxStrain = strain;
      maxIndex = i;
    }
  }
  return maxIndex;
}

void networkV4::propogator::breakMostStrained(network& _network, bondType _type)
{
  bonds& bonds = _network.getBonds();
  size_t maxIndex = getMostStrained(_network, m_breakType);
  bond& b = bonds.get(maxIndex);
  b.connected() = false;
  m_dataOut->writeBondData(genBondData(_network, b));
}

auto networkV4::propogator::findSingleBreak(network& _network) -> bool
{
  network testNetwork;
  bool converged = false;
  double distToBreakA, distToBreakB;
  size_t mostStrainedIndex;
  bond mostStrainedBond;

  double a = getStrain(_network);
  double b = getStrain(_network) + m_maxStep;

  testNetwork = _network;
  evalStrain(testNetwork, a);
  mostStrainedIndex = getMostStrained(testNetwork, m_breakType);
  mostStrainedBond = testNetwork.getBonds().get(mostStrainedIndex);
  distToBreakA =
      mostStrainedBond.lambda() - testNetwork.bondStrain(mostStrainedBond);

  testNetwork = _network;
  evalStrain(testNetwork, b);
  mostStrainedIndex = getMostStrained(testNetwork, m_breakType);
  mostStrainedBond = testNetwork.getBonds().get(mostStrainedIndex);
  distToBreakB =
      mostStrainedBond.lambda() - testNetwork.bondStrain(mostStrainedBond);

  if (distToBreakB > 0.0) {
    _network = testNetwork;
    return false;
  }

  roots::ITP solver(a, b, m_tol);
  for (size_t iters = 0; iters < solver.nMax(); ++iters) {
    double xITP = solver.guessRoot(a, b, distToBreakA, distToBreakB);

    double distToBreakITP;
    testNetwork = _network;
    evalStrain(testNetwork, xITP);
    mostStrainedIndex = getMostStrained(testNetwork, m_breakType);
    mostStrainedBond = testNetwork.getBonds().get(mostStrainedIndex);
    distToBreakITP =
        mostStrainedBond.lambda() - testNetwork.bondStrain(mostStrainedBond);
    if (distToBreakITP >= 0.) {
      b = xITP;
      distToBreakB = distToBreakITP;
    } else {
      a = xITP;
      distToBreakB = distToBreakITP;
    }
    if (std::abs(b - a) < 2 * m_tol) {
      break;
    }
  }
  _network = testNetwork;
  return true;
}

auto networkV4::propogator::genTimeData(const network& _network,
                                        const std::string& _reason,
                                        size_t _breakCount)
    -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(10 + 5 * types.size());
  data.emplace_back(_reason);
  data.emplace_back(m_strainCount);
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

auto networkV4::propogator::genBondData(
    const network& _network, const bond& _bond) -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(23 + 5 * types.size());
  data.emplace_back(m_strainCount);
  data.emplace_back(enumString::bondType2Str.at(_bond.type()));
  data.emplace_back(_bond.mu());
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