#include "Minimization.hpp"

networkV4::SteepestDescent::SteepestDescent() {}

networkV4::SteepestDescent::SteepestDescent(double _tol)
    : m_tol(_tol)
{
}

networkV4::SteepestDescent::~SteepestDescent() {}

auto networkV4::SteepestDescent::integrate(network& _network) -> void
{
  _network.computeForces();
  for (size_t iter = 0; iter < m_maxIter; ++iter) {
    m_dt = lineSearch(_network);
    Integrationtools::overdampedMove(_network, m_dt);
    _network.computeForces();
    double error = tools::norm(_network.getNodes().forces());
    if (error < m_tol) {
      break;
    }
  }
}

auto networkV4::SteepestDescent::getDt() const -> double
{
  return m_dt;
}

auto networkV4::SteepestDescent::eGradient(network& _network,
                                           const double _step) const -> double
{
  /// const double norm = tools::norm(_network.getNodes().forces());
  Integrationtools::overdampedMove(_network, -_step);
  const double energy = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, 2.0 * _step);
  const double newEnergy = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, -_step);
  return (newEnergy - energy) / (2.0 * _step);
}

auto networkV4::SteepestDescent::lineSearch(network& _network) const -> double
{
  double a = m_dtMin;
  double b = m_dtMax;

  double initialEnergy = _network.computeEnergy();
  // overdampedMove(_network, a);
  double fa = eGradient(_network);
  // double Ea = _network.computeEnergy();
  // overdampedMove(_network, -a);
  Integrationtools::overdampedMove(_network, b);
  double fb = eGradient(_network);
  double Eb = _network.computeEnergy();
  Integrationtools::overdampedMove(_network, -b);

  roots::ITP solver(a, b, m_rootTol);
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

    if (std::abs(b - a) < 2 * m_rootTol) {
      break;
    }
  }
  return (a + b) * 0.5;
}

networkV4::FireMinimizer::FireMinimizer() {}

networkV4::FireMinimizer::FireMinimizer(double _tol)
    : m_tol(_tol)
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
  size_t Npos = 0;
  size_t Nneg = 0;
  double alpha = m_alpha0;

  double scale1, scale2, abc, vdotf, vdotv, fdotf;

  nodes& networkNodes = _network.getNodes();
  auto& vels = networkNodes.velocities();
  auto& forces = networkNodes.forces();
  _network.computeForces();

  double error = tools::norm(forces);
  if (error < m_tol) {
    return;
  }

  networkNodes.clearVelocities();
  Integrationtools::updateVelocities(_network, m_dt);

  for (size_t iter = 0; iter < m_maxIter; ++iter) {
    vdotf = tools::xdoty(vels, forces);
    if (vdotf > 0.0) {
      Npos++;
      Nneg = 0;
      vdotv =
          tools::xdoty(networkNodes.velocities(), networkNodes.velocities());
      fdotf = tools::xdoty(networkNodes.forces(), networkNodes.forces());

      scale1 = 1.0 - alpha;
      scale2 = fdotf <= 1e-20 ? 0.0 : (alpha * std::sqrt(vdotv / fdotf));

      if (Npos > m_Ndelay) {
        m_dt = std::min(m_dt * m_finc, m_dtMax);
        alpha *= m_falpha;
      }
    } else {
      Nneg++;
      Npos = 0;

      if (Nneg > m_Nnegmax) {
        break;
      }
      if (iter > m_Ndelay) {
        m_dt = std::max(m_dt * m_fdec, m_dtMin);
        alpha = m_alpha0;
      }
      Integrationtools::move(_network, -0.5 * m_dt);
      networkNodes.clearVelocities();
    }
    Integrationtools::updateVelocities(_network, m_dt);
    double maxAbsVel = tools::maxAbsComponent(vels);
    m_dt = std::max(std::min(m_dt, m_dmax / maxAbsVel), m_dtMin);

    if (vdotf > 0.0) {
      const double vmax = m_dmax / m_dt;
      for (size_t i = 0; i < vels.size(); ++i) {
        vels[i] = (vels[i] * scale1) + (forces[i] * scale2);
      }
    }
    Integrationtools::move(_network, m_dt);
    _network.computeForces();

    // double test = tools::maxLength(forces);
    // double test2 = tools::maxAbsComponent(forces);
    error = tools::norm(forces);
    if (error < m_tol) {
      break;
    }
  }
}

networkV4::OverdampedAdaptiveMinimizer::OverdampedAdaptiveMinimizer() {}

networkV4::OverdampedAdaptiveMinimizer::OverdampedAdaptiveMinimizer(double _esp,
                                                                    double _tol)
    : m_esp(_esp)
    , m_tol(_tol)
{
}

networkV4::OverdampedAdaptiveMinimizer::~OverdampedAdaptiveMinimizer() {}

void networkV4::OverdampedAdaptiveMinimizer::integrate(network& _network)
{
  double nextDt = m_dt;
  bool error = false;

  nodes& networkNodes = _network.getNodes();
  for (size_t iter = 0; iter < m_maxOuterIter; ++iter) {
    nextDt = innerIteration(_network, nextDt);
    double error = tools::norm(networkNodes.forces());
    if (error < m_tol) {
      break;
    }
    m_dt = nextDt;
  }
}

auto networkV4::OverdampedAdaptiveMinimizer::innerIteration(
    network& _network, double _dt) -> double
{
  m_dt = _dt;
  double q = m_qMin;
  bool error = false;
  size_t iter = 0;
  nodes& networkNodes = _network.getNodes();
  double startEnergy = _network.getEnergy();

  const double maxAbsVel = tools::maxAbsComponent(networkNodes.forces());
  m_dt = std::max(std::min(m_dt, m_dmax / maxAbsVel), m_dtMin);

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

    const double forceError =
        Integrationtools::vectorDiffNorm(networkNodes.forces(),
                                         m_frnbar,
                                         networkNodes.fixed(),
                                         true);  // |f(rbar_{n+1}) - f(r_n)|
    q = std::clamp(std::pow(m_esp / forceError, 2), m_qMin, m_qMax);
    const double energy = _network.getEnergy();

    if (q > 1.0 && energy < startEnergy)
      break;

    error = (m_dt == m_dtMin || std::isnan(q));
    if (error)
      break;

    if (energy > startEnergy)
        q = std::min(q, m_fdec);
    q = std::clamp(q, m_qMin, m_qMax);
    
    Integrationtools::overdampedMove(_network, -m_dt);  // r_{n+1} = r_n
    Integrationtools::copyVector(m_frn, networkNodes.forces());
    m_dt = std::clamp(m_dt * q, m_dtMin, m_dtMax);
  }
  if (error) {
    throw std::runtime_error("Adaptive Euler Heun failed to converge");
  }
  return std::clamp(m_dt * q, m_dtMin, m_dtMax);
}

auto networkV4::OverdampedAdaptiveMinimizer::getDt() const -> double
{
  return m_dt;
}