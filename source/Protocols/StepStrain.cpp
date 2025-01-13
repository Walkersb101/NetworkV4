#include "StepStrain.hpp"

#include "Misc/Tools.hpp"
#include "Integration/OverDamped.hpp"

networkV4::stepStrain::stepStrain() {}

networkV4::stepStrain::stepStrain(double _maxStrain,
                                  StrainType _strainType,
                                  std::unique_ptr<BreakTypes>& _breakType,
                                  double _maxTime,
                                  double _esp,
                                  double _timeScale,
                                  double _stressScale)
    : m_maxStrain(_maxStrain)
    , m_maxTime(_maxTime)
    , m_esp(_esp)
    , m_t(0.0)
    , protocol(_strainType, _breakType)
    , m_timeScale(_timeScale)
    , m_stressScale(_stressScale)
    , m_logCount(2)
{
}

networkV4::stepStrain::~stepStrain() {}

void networkV4::stepStrain::run(network& _network)
{
  OverdampedAdaptiveEulerHeun integrator(config::integrators::default_dt,
                                         m_esp);

  m_networkOut->save(_network, 0, 0.0, "Initial");
  m_dataOut->writeTimeData(genTimeData(_network, "Initial"));

  strain(_network, m_maxStrain);
  _network.computeForces();
  m_networkOut->save(_network, 1, 0.0, "Affine");
  m_dataOut->writeTimeData(genTimeData(_network, "Affine"));

  const double initialStress = globalStressAlongAxis(_network);
  m_stressLogPar = std::abs(initialStress) * m_stressScale;
  m_timeLogPar = integrator.getDt() * m_timeScale;

  while (m_t < m_maxTime) {
    integrator.integrate(_network);
    m_t += integrator.getDt();
    const auto broken = m_breakProtocol->Break(_network, integrator);
    if (broken.size() > 0) {
        _network.computeForces();
    }

    for (const auto b : broken) {
      m_dataOut->writeBondData(genBondData(_network, *b, m_breakProtocol));
    }
    checkLog(_network);


    if (checkExit(_network)) {
      break;
    }
  }
  m_networkOut->save(_network, m_logCount, m_t, "Final");
  m_dataOut->writeTimeData(genTimeData(_network, "Final"));
}

void networkV4::stepStrain::initIO(const network& _network,
                                   std::unique_ptr<dataOut>& _dataOut,
                                   std::unique_ptr<networkOut>& _networkOut)
{
  auto types = tools::uniqueBondTypes(_network.getBonds());

  std::vector<std::string> dataHeader = {
      "Reason",
      "Time",
      "Domainx",
      "Domainy",
      enumString::strainType2Str.at(m_strainType) + "Strain",
      "RMSVel",
      "MaxVel"};

  std::vector<std::string> bondHeader = {
      "Time",
      "Type",
      "mu",
      "E",
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
      "MaxForce",
      "RMSVel",
      "MaxVel"};

  tensorData::addBondCountByTypeHeader(dataHeader, types);
  tensorData::addBondCountByTypeHeader(bondHeader, types);

  tensorData::addStressByTypeHeader(bondHeader, types);
  tensorData::addStressByTypeHeader(dataHeader, types);

  _dataOut->initFiles(dataHeader, bondHeader);
  m_dataOut = _dataOut.get();
  m_networkOut = _networkOut.get();
}

auto networkV4::stepStrain::genTimeData(const network& _network,
                                        const std::string& _reason)
    -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(11 + 5 * types.size());
  data.emplace_back(_reason);
  data.emplace_back(m_t);
  data.emplace_back(_network.getDomain()[0]);
  data.emplace_back(_network.getDomain()[1]);
  data.emplace_back(getStrain(_network));
  data.emplace_back(tools::norm(_network.getNodes().velocities()));
  data.emplace_back(tools::maxLength(_network.getNodes().velocities()));

  data.emplace_back(_network.getBonds().connectedCount());
  for (const auto& type : types) {
    data.emplace_back(_network.getBonds().connectedCount(type));
  }

  data.insert(
      data.end(),
      Utils::Math::flatten(globalStress));
  for (const auto& type : types) {
    tensor2 stress = stresses[type];
    data.insert(data.end(), Utils::Math::flatten(globalStress));
  }
  return data;
}

auto networkV4::stepStrain::genBondData(const network& _network,
                                        const bond& _bond,
                                        BreakTypes* _breakProtocol)
    -> std::vector<writeableTypes>
{
  auto types = tools::uniqueBondTypes(_network.getBonds());
  auto stresses = _network.getStresses();
  auto globalStress = _network.getGlobalStress();

  std::vector<writeableTypes> data;
  data.reserve(24 + 5 * types.size());
  data.emplace_back(m_t);
  data.emplace_back(enumString::bondType2Str.at(_bond.type()));
  data.emplace_back(_bond.mu());
  data.emplace_back(_bond.lambda());
  data.emplace_back(_bond.naturalLength());
  data.emplace_back(_network.bondStrain(_bond));
  data.emplace_back(_bond.src());
  data.emplace_back(_bond.dst());
  data.emplace_back(_network.getNodes().position(_bond.src())[0]);
  data.emplace_back(_network.getNodes().position(_bond.src())[1]);
  data.emplace_back(_network.getNodes().position(_bond.dst())[0]);
  data.emplace_back(_network.getNodes().position(_bond.dst())[1]);
  data.emplace_back(_network.getDomain()[0]);
  data.emplace_back(_network.getDomain()[1]);
  data.emplace_back(getStrain(_network));
  data.emplace_back(tools::norm(_network.getNodes().forces()));
  data.emplace_back(tools::maxLength(_network.getNodes().forces()));
  data.emplace_back(tools::norm(_network.getNodes().velocities()));
  data.emplace_back(tools::maxLength(_network.getNodes().velocities()));

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

double networkV4::stepStrain::globalStressAlongAxis(
    const network& _network) const
{
  tensor2 stress = _network.getGlobalStress();
  switch (m_strainType) {
    case StrainType::Shear:
      return stress.xy;
    case StrainType::Elongation:
      return stress.yy;
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

void networkV4::stepStrain::checkLog(const network& _network)
{
  std::string reason;
  double stress = globalStressAlongAxis(_network);
  if (m_t > m_timeLogPar) {
    reason = "Time";
    m_timeLogPar = std::max(m_timeLogPar * m_timeScale, m_t);
  } else if (std::abs(stress) < m_stressLogPar) {
    reason = "Stress";
    m_stressLogPar = std::min(m_stressLogPar * m_stressScale, stress);
  } else {
    return;
  }
  m_dataOut->writeTimeData(genTimeData(_network, reason));
  m_networkOut->save(
      _network, m_logCount++, m_t, reason + "-" + std::to_string(m_logCount));
}

auto networkV4::stepStrain::checkExit(const network& _network) -> bool
{
  double stress = globalStressAlongAxis(_network);
  return std::abs(stress) <= 1e-10;
}