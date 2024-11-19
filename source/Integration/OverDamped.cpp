#include "OverDamped.hpp"

#include "Misc/Tools.hpp"

networkV4::overdampedEuler::overdampedEuler()
    : m_dt(config::integrators::default_dt)
{
}

networkV4::overdampedEuler::overdampedEuler(double _dt)
    : m_dt(_dt)
{
}

networkV4::overdampedEuler::~overdampedEuler() {}

void networkV4::overdampedEuler::integrate(network& _network)
{
  _network.computeForces();
  Integrationtools::overdampedMove(_network, m_dt);
}

auto networkV4::overdampedEuler::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedEulerHeun::OverdampedEulerHeun()
    : m_dt(config::integrators::default_dt)
{
}

networkV4::OverdampedEulerHeun::OverdampedEulerHeun(double _dt)
    : m_dt(_dt)
{
}

networkV4::OverdampedEulerHeun::~OverdampedEulerHeun() {}

void networkV4::OverdampedEulerHeun::integrate(network& _network)
{
  _network.computeForces();

  m_tempForces = _network.getNodes().forces();
  m_tempStresses = _network.getStresses();

  Integrationtools::overdampedMove(_network, m_dt);
  _network.computeForces();
  Integrationtools::heunAverage(_network, m_tempForces, m_dt);

  _network.getStresses() =
      tools::averageMaps(_network.getStresses(), m_tempStresses);
}

auto networkV4::OverdampedEulerHeun::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedAdaptiveEulerHeun::OverdampedAdaptiveEulerHeun()
    : m_dt(config::integrators::default_dt)
    , m_nextdt(config::integrators::default_dt)
    , m_esp(config::integrators::adaptiveIntegrator::esp)
{
}

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
  double q;

  std::size_t iters = 0;

  std::size_t maxIter = config::integrators::adaptiveIntegrator::maxIter;
  double qMin = config::integrators::adaptiveHeun::qMin;
  double qMax = config::integrators::adaptiveHeun::qMax;
  double dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double dtMax = config::integrators::adaptiveIntegrator::dtMax;

  nodes& networkNodes = _network.getNodes();

  _network.computeForces();
  m_tempForces = networkNodes.forces();
  m_tempPositions = networkNodes.positions();
  m_tempStresses = _network.getStresses();

  while (iters < maxIter) {
    Integrationtools::overdampedMove(_network, m_dt);
    _network.computeForces();

    const double posError = forceErrorNorm(networkNodes) * m_dt;
    q = std::clamp(std::pow(m_esp / posError, 2), qMin, qMax);
    if (q > 1.0) {
      break;
    }
    _network.getNodes().positions() = m_tempPositions;
    m_dt *= q;
    ++iters;
  }
  if (iters == maxIter) {
    std::runtime_error("Max iterations reached");
  }

  Integrationtools::heunAverage(_network, m_tempForces, m_dt);

  _network.getStresses() =
      tools::averageMaps(_network.getStresses(), m_tempStresses);

  m_nextdt = std::clamp(m_dt * q, dtMin, dtMax);
}

auto networkV4::OverdampedAdaptiveEulerHeun::getDt() const -> double
{
  return m_dt;
}

auto networkV4::OverdampedAdaptiveEulerHeun::forceErrorNorm(nodes& _nodes)
    -> double
{
  double sum = 0.0;
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    const double force = (m_tempForces[i] - _nodes.force(i)).lengthSquared();
    sum += force;
  }
  return sqrt(sum);
}