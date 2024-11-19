#include "Minimization.hpp"

networkV4::SteepestDescent::SteepestDescent(double _dt)
    : m_dt(_dt)
    , m_nextdt(_dt)
    , m_prevEnergy(std::numeric_limits<double>::quiet_NaN())
{
}

networkV4::SteepestDescent::SteepestDescent()
    : SteepestDescent(config::integrators::default_dt)
{
}

networkV4::SteepestDescent::~SteepestDescent() {}

auto networkV4::SteepestDescent::integrate(network& _network) -> void
{
  m_dt = m_nextdt;
  double energy;

  std::size_t maxIter = config::integrators::adaptiveIntegrator::maxIter;
  double qMin = config::integrators::adaptiveHeun::qMin;
  double qMax = config::integrators::adaptiveHeun::qMax;
  double dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double dtMax = config::integrators::adaptiveIntegrator::dtMax;

  double gamma = 1e-4;

  std::size_t iters = 0;

  nodes& networkNodes = _network.getNodes();
  _network.computeForces();
  if (std::isnan(m_prevEnergy)) {
    m_prevEnergy = _network.getEnergy();
  }
  // const double maxLen = tools::maxLength(networkNodes.forces());
  // const double FdotF =
  //     tools::xdoty(networkNodes.forces(), networkNodes.forces());

  m_dt = lineSearch(_network);
  // m_dt = std::clamp(m_dt, dtMin, dtMax);
  Integrationtools::overdampedMove(_network, m_dt);
  energy = _network.computeEnergy();
}

auto networkV4::SteepestDescent::getDt() const -> double
{
  return m_dt;
}

auto networkV4::SteepestDescent::eGradient(network& _network,
                                           const double _step) const -> double
{
  ///const double norm = tools::norm(_network.getNodes().forces());
  Integrationtools::overdampedMove(_network, - _step);
  const double energy = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, 2.0 * _step);
  const double newEnergy = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, -_step);
  return (newEnergy - energy) / ( 2.0 * _step);
}

auto networkV4::SteepestDescent::lineSearch(network& _network) const -> double
{
  double a = config::integrators::adaptiveIntegrator::dtMin;
  double b = config::integrators::adaptiveIntegrator::dtMax;
  double tol = 1e-8;

  double initialEnergy = _network.computeEnergy();
  // overdampedMove(_network, a);
  double fa = eGradient(_network);
  // double Ea = _network.computeEnergy();
  // overdampedMove(_network, -a);
  Integrationtools::overdampedMove(_network, b);
  double fb = eGradient(_network);
  double Eb = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, -b);

  roots::ITP solver(a, b, tol);
  if (fa * fb > 0.0) {
    if (Eb < initialEnergy) {
      return b;
    } else {
      return a;
    }
  }

  for (std::size_t iters = 0; iters < solver.nMax(); ++iters) {
    double xITP = solver.guessRoot(a, b, fa, fb);

    Integrationtools::overdampedMove(_network, xITP);
    double fITP = eGradient(_network);
    Integrationtools::overdampedMove(_network, -xITP);
    if (fITP >= 0.) {
      b = xITP;
      fb = fITP;
    } else {
      a = xITP;
      fa = fITP;
    }

    if (std::abs(b - a) < 2 * tol) {
      break;
    }
  }
  return (a + b) * 0.5;
}

networkV4::FireMinimizer::FireMinimizer()
    : m_dt(config::integrators::default_dt)
    , m_tol(config::integrators::miminizer::tol)
{
}

networkV4::FireMinimizer::FireMinimizer(double _tol)
    : m_dt(config::integrators::default_dt)
    , m_tol(_tol)
{
}

networkV4::FireMinimizer::FireMinimizer(double _dt, double _tol)
    : m_dt(_dt)
    , m_tol(_tol)
{
}

networkV4::FireMinimizer::~FireMinimizer() {}

void networkV4::FireMinimizer::integrate(network& _network)
{
  size_t maxIter = config::integrators::miminizer::maxIter;
  double dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double dtMax = config::integrators::adaptiveIntegrator::dtMax;

  double alpha0 = config::integrators::fire2::alpha0;
  double finc = config::integrators::fire2::finc;
  double fdec = config::integrators::fire2::fdec;
  double falpha = config::integrators::fire2::falpha;
  size_t Ndelay = config::integrators::fire2::Ndelay;
  size_t Nnegmax = config::integrators::fire2::Nnegmax;
  double dmax = config::integrators::fire2::dmax;

  size_t Npos = 0;
  size_t Nneg = 0;
  double alpha = alpha0;

  double scale1, scale2, abc, vdotf, vdotv, fdotf;

  nodes& networkNodes = _network.getNodes();
  _network.computeForces();

  double error = tools::norm(networkNodes.forces());
  if (error < m_tol) {
    return;
  }

  networkNodes.clearVelocities();
  Integrationtools::updateVelocities(_network, m_dt);

  for (size_t iter = 0; iter < maxIter; ++iter) {
    vdotf = tools::xdoty(networkNodes.velocities(), networkNodes.forces());
    if (vdotf > 0.0) {
      Npos++;
      Nneg = 0;
      vdotv =
          tools::xdoty(networkNodes.velocities(), networkNodes.velocities());
      fdotf = tools::xdoty(networkNodes.forces(), networkNodes.forces());

      scale1 = 1.0 - alpha;
      scale2 = fdotf <= 1e-20 ? 0.0 : (alpha * std::sqrt(vdotv / fdotf));

      if (Npos > Ndelay) {
        m_dt = std::min(m_dt * finc, dtMax);
        alpha *= falpha;
      }
    } else {
      Nneg++;
      Npos = 0;

      if (Nneg > Nnegmax) {
        break;
      }
      if (iter > Ndelay) {
        m_dt = std::max(m_dt * fdec, dtMin);
        alpha = alpha0;
      }
      Integrationtools::move(_network, -0.5 * m_dt);
      networkNodes.clearVelocities();
      _network.computeForces();
    }
    Integrationtools::updateVelocities(_network, m_dt);
    double maxAbsVel = tools::maxAbsComponent(networkNodes.velocities());
    m_dt = std::min(m_dt, dmax / maxAbsVel);

    if (vdotf > 0.0) {
      auto& vels = networkNodes.velocities();
      auto& forces = networkNodes.forces();
      const double vmax = dmax / m_dt;
      for (size_t i = 0; i < vels.size(); ++i) {
        vels[i] = (vels[i] * scale1) + (forces[i] * scale2);
      }
    }
    Integrationtools::move(_network, m_dt);
    _network.computeForces();

    double test = tools::maxLength(networkNodes.forces());
    double test2 = tools::maxAbsComponent(networkNodes.forces());
    error = tools::norm(networkNodes.forces());
    if (error < m_tol) {
      break;
    }
  }
}

networkV4::OverdampedAdaptiveMinimizer::OverdampedAdaptiveMinimizer()
    : m_integrator()
{
}

networkV4::OverdampedAdaptiveMinimizer::OverdampedAdaptiveMinimizer(double _dt,
                                                                    double _esp,
                                                                    double _tol)
    : m_integrator(_dt, _esp)
    , m_tol(_tol)
{
}

networkV4::OverdampedAdaptiveMinimizer::~OverdampedAdaptiveMinimizer() {}

void networkV4::OverdampedAdaptiveMinimizer::integrate(network& _network)
{
  size_t maxIter = config::integrators::miminizer::maxIter;

  size_t iters = 0;
  while (iters < maxIter) {
    m_integrator.integrate(_network);
    double error = tools::norm(_network.getNodes().forces());
    if (error < m_tol) {
      break;
    }
    ++iters;
  };
}
