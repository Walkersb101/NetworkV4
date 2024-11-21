#include "OverDamped.hpp"

#include "Misc/Tools.hpp"

networkV4::overdampedEuler::overdampedEuler() {}

networkV4::overdampedEuler::overdampedEuler(double _dt)
    : m_dt(_dt)
{
}

networkV4::overdampedEuler::~overdampedEuler() {}

void networkV4::overdampedEuler::integrate(network& _network)
{
  Integrationtools::overdampedMove(_network, m_dt);
  _network.computeForces();
}

auto networkV4::overdampedEuler::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedEulerHeun::OverdampedEulerHeun() {}

networkV4::OverdampedEulerHeun::OverdampedEulerHeun(double _dt)
    : m_dt(_dt)
{
}

networkV4::OverdampedEulerHeun::~OverdampedEulerHeun() {}

void networkV4::OverdampedEulerHeun::integrate(network& _network)
{
  Integrationtools::copyVector(_network.getNodes().forces(), m_frn);

  Integrationtools::overdampedMove(_network,
                                   m_dt);  // rbar_{n+1} = r_n + f(r_n) * dt
  _network.computeForces();
  Integrationtools::heunAverage(
      _network,
      m_frn,
      m_dt);  // r_{n+1} = rn + 0.5 * (f(r_n) + f(rbar_{n+1})) * dt
              // = rbar_{n+1} + 0.5 * (f(rbar_{n+1}) - f(r_n)) * dt
  _network.computeForces();
}

auto networkV4::OverdampedEulerHeun::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedAdaptiveEulerHeun::OverdampedAdaptiveEulerHeun() {}

networkV4::OverdampedAdaptiveEulerHeun::OverdampedAdaptiveEulerHeun(double _dt,
                                                                    double _esp)
    : m_dt(_dt)
    , m_nextdt(_dt)
    , m_esp(_esp)
{
}

networkV4::OverdampedAdaptiveEulerHeun::~OverdampedAdaptiveEulerHeun() {}

void networkV4::OverdampedAdaptiveEulerHeun::integrate(network& _network)
{
  m_dt = m_nextdt;
  double q = m_qMin;
  bool error = false;
  size_t iter = 0;
  nodes& networkNodes = _network.getNodes();

  Integrationtools::copyVector(networkNodes.forces(), m_frn);
  while (iter++ < m_maxInnerIter) {
    Integrationtools::overdampedMove(_network,
                                     m_dt);  // rbar_{n+1} = r_n + f(r_n) * dt
    _network.computeForces();
    Integrationtools::copyVector(networkNodes.forces(),
                                 m_frnbar);  // f(rbar_{n+1})
    Integrationtools::heunAverage(
        _network,
        m_frn,  // r_{n+1} = rn + 0.5 * (f(r_n) + f(rbar_{n+1})) * dt
        m_dt);  // = rbar_{n+1} + 0.5 * (f(rbar_{n+1}) - f(r_n)) * dt
    _network.computeForces();

    const double forceError = Integrationtools::vectorDiffNorm(
        networkNodes.forces(), m_frnbar,
        networkNodes.fixed(), true);  // |f(rbar_{n+1}) - f(r_n)|
    const double q = std::clamp(std::pow(m_esp / forceError, 2), m_qMin, m_qMax);
    
    if (q > 1.0)
      break;
    error = (m_dt == m_dtMin || std::isnan(q));
    if (error)
      break;

    Integrationtools::overdampedMove(_network, -m_dt);  // r_{n+1} = r_n
    Integrationtools::copyVector(m_frn, networkNodes.forces());
    m_dt = std::clamp(m_dt * q, m_dtMin, m_dtMax);
  }
  if (error) {
    throw std::runtime_error("Adaptive Euler Heun failed to converge");
  }
  m_nextdt = std::clamp(m_dt * q, m_dtMin, m_dtMax);
}

auto networkV4::OverdampedAdaptiveEulerHeun::getDt() const -> double
{
  return m_dt;
}