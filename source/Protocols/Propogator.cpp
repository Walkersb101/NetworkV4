#include "Propogator.hpp"

#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"
#include "Integration/Minimization.hpp"

networkV4::propogator::propogator() {}

networkV4::propogator::propogator(std::vector<double> _strains,
                                  StrainType _strainType)
    : propogator(_strains,
                 _strainType,
                 bondType::any,
                 config::integrators::adaptiveIntegrator::esp,
                 config::integrators::miminizer::tol,
                 0.0)
{
}

networkV4::propogator::propogator(const std::vector<double>& _strains,
                                  StrainType _strainType,
                                  bondType _breakType,
                                  double _esp,
                                  double _tol,
                                  double _maxStep)
    : m_strains(_strains)
    , m_breakType(_breakType)
    , m_esp(_esp)
    , m_tol(_tol)
    , m_maxStep(_maxStep)
    , protocol(_strainType)
    , m_strainCount(0)
{
}

networkV4::propogator::~propogator() {}

void networkV4::propogator::run(network& _network)
{
  relax(_network);
  m_dataOut->writeTimeData(genTimeData(_network, "Initial", 0));
  m_networkOut->save(_network, 0, 0.0, "Initial");

  network InitalNetwork = _network;
  for (const auto& strain : m_strains) {
    m_strainCount++;
    evalStrain(_network, strain);
    if (m_breakType == bondType::any)
      breakMostStrained(_network);
    else
      breakMostStrained(_network, m_breakType);
    m_dataOut->writeTimeData(genTimeData(_network, "Start", 1));
    m_networkOut->save(_network, m_strainCount, 0.0, "Start");
    relax(_network);
    m_dataOut->writeTimeData(genTimeData(_network, "End", 1));
    m_networkOut->save(_network, m_strainCount, 1.0, "End");
    _network = _network;
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

void networkV4::propogator::evalStrain(network& _network, double _strain)
{
  if (m_maxStep == 0.0) {
    strain(_network, _strain - getStrain(_network));
    _network.computeForces();
    relax(_network);
    return;
  }
  while (getStrain(_network) < _strain - 1e-15) {
    const double strainStep =
        std::min(_strain - getStrain(_network), m_maxStep);
    strain(_network, strainStep);
    _network.computeForces();
    relax(_network);
  }
}

void networkV4::propogator::breakMostStrained(network& _network)
{
  bonds& bonds = _network.getBonds();
  std::vector<double> strains;
  strains.reserve(bonds.size());
  for (const auto& b : bonds) {
    if (!b.connected()) {
      strains.push_back(0.0);
      continue;
    }
    const vec2d bDist = _network.minDist(_network.getNodes().position(b.src()),
                                         _network.getNodes().position(b.dst()));
    const double strain =
        (bDist.length() - b.naturalLength()) / b.naturalLength();
    strains.push_back(strain);
  }
  auto maxElementIt = std::max_element(strains.begin(), strains.end());
  size_t maxIndex = std::distance(strains.begin(), maxElementIt);
  bond& b = bonds.get(maxIndex);
  b.connected() = false;
  m_dataOut->writeBondData(genBondData(_network, maxIndex));
}

void networkV4::propogator::breakMostStrained(network& _network, bondType _type)
{
  bonds& bonds = _network.getBonds();
  std::vector<double> strains;
  strains.reserve(bonds.size());
  for (const auto& b : bonds) {
    if (!b.connected() || b.type() != _type) {
      strains.push_back(0.0);
      continue;
    }
    const vec2d bDist = _network.minDist(_network.getNodes().position(b.src()),
                                         _network.getNodes().position(b.dst()));
    const double strain =
        (bDist.length() - b.naturalLength()) / b.naturalLength();
    strains.push_back(strain);
  }
  auto maxElementIt = std::max_element(strains.begin(), strains.end());
  size_t maxIndex = std::distance(strains.begin(), maxElementIt);
  bond& b = bonds.get(maxIndex);
  b.connected() = false;
  m_dataOut->writeBondData(genBondData(_network, maxIndex));
}

void networkV4::propogator::relax(network& _network)
{
  FireMinimizer minimizer(m_tol);
  minimizer.integrate(_network);
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
    const network& _network, size_t _bondIndex) -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();
  auto& b = _network.getBonds().get(_bondIndex);

  std::vector<writeableTypes> data;
  data.reserve(23 + 5 * types.size());
  data.emplace_back(m_strainCount);
  data.emplace_back(enumString::bondType2Str.at(b.type()));
  data.emplace_back(b.mu());
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