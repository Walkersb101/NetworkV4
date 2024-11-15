#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "Intergrator.hpp"

#include "Misc/Config.hpp"
#include "Misc/Tools.hpp"
#include "Core/Network.hpp"
#include "Core/Nodes.hpp"
#include "Core/Roots.hpp"
#include "Core/Tensor2.hpp"
#include "Core/Vec2.hpp"

void overdampedMove(networkV4::network& _network, double _dt)
{
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.position(i) += _nodes.force(i) * _dt;
    _network.wrapPosition(_nodes.position(i));
  }
}

void move(networkV4::network& _network, double _dt)
{
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.position(i) += _nodes.velocity(i) * _dt;
    _network.wrapPosition(_nodes.position(i));
  }
}

void updateVelocities(networkV4::network& _network, double _alpha)
{
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.velocity(i) += _nodes.force(i) * _alpha / _nodes.mass(i);
  }
}

void heunAverage(networkV4::network& _network,
                 const std::vector<vec2d>& _tempForces,
                 double _dt)
{
  networkV4::nodes& networkNodes = _network.getNodes();
  for (std::size_t i = 0; i < networkNodes.size(); ++i) {
    if (networkNodes.fixed(i)) {
      continue;
    }
    networkNodes.position(i) +=
        (networkNodes.force(i) - _tempForces[i]) * _dt * 0.5;
    _network.wrapPosition(networkNodes.position(i));
    networkNodes.force(i) = (networkNodes.force(i) + _tempForces[i]) * 0.5;
  }
}

networkV4::intergrator::intergrator() {}

networkV4::intergrator::~intergrator() {}

void networkV4::intergrator::integrate(network& _network) {}

auto networkV4::intergrator::getDt() const -> double
{
  return 0.0;
}

networkV4::overdampedEuler::overdampedEuler()
    : m_dt(config::intergrators::default_dt)
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
  overdampedMove(_network, m_dt);
}

auto networkV4::overdampedEuler::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedEulerHeun::OverdampedEulerHeun()
    : m_dt(config::intergrators::default_dt)
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

  overdampedMove(_network, m_dt);
  _network.computeForces();
  heunAverage(_network, m_tempForces, m_dt);

  _network.getStresses() =
      tools::averageMaps(_network.getStresses(), m_tempStresses);
}

auto networkV4::OverdampedEulerHeun::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedAdaptiveEulerHeun::OverdampedAdaptiveEulerHeun()
    : m_dt(config::intergrators::default_dt)
    , m_nextdt(config::intergrators::default_dt)
    , m_esp(config::intergrators::adaptiveIntergrator::esp)
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

  std::size_t maxIter = config::intergrators::adaptiveIntergrator::maxIter;
  double qMin = config::intergrators::adaptiveHeun::qMin;
  double qMax = config::intergrators::adaptiveHeun::qMax;
  double dtMin = config::intergrators::adaptiveIntergrator::dtMin;
  double dtMax = config::intergrators::adaptiveIntergrator::dtMax;

  nodes& networkNodes = _network.getNodes();

  _network.computeForces();
  m_tempForces = networkNodes.forces();
  m_tempPositions = networkNodes.positions();
  m_tempStresses = _network.getStresses();

  while (iters < maxIter) {
    overdampedMove(_network, m_dt);
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

  heunAverage(_network, m_tempForces, m_dt);

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

networkV4::SteepestDescent::SteepestDescent(double _dt)
    : m_dt(_dt)
    , m_nextdt(_dt)
    , m_prevEnergy(std::numeric_limits<double>::quiet_NaN())
{
}

networkV4::SteepestDescent::SteepestDescent()
    : SteepestDescent(config::intergrators::default_dt)
{
}

networkV4::SteepestDescent::~SteepestDescent() {}

auto networkV4::SteepestDescent::integrate(network& _network) -> void
{
  m_dt = m_nextdt;
  double energy;

  std::size_t maxIter = config::intergrators::adaptiveIntergrator::maxIter;
  double qMin = config::intergrators::adaptiveHeun::qMin;
  double qMax = config::intergrators::adaptiveHeun::qMax;
  double dtMin = config::intergrators::adaptiveIntergrator::dtMin;
  double dtMax = config::intergrators::adaptiveIntergrator::dtMax;

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
  overdampedMove(_network, m_dt);
  energy = _network.computeEnergy();

  // while (energy > m_prevEnergy - gamma * m_dt * FdotF && iters < maxIter
  //        && m_dt > dtMin)
  //{
  //   overdampedMove(_network, -m_dt);
  //   m_dt = ((m_dt * m_dt) * FdotF)
  //       / (2.0 * (energy + m_dt * FdotF - m_prevEnergy));
  //   m_dt = std::clamp(m_dt, dtMin, dtMax);
  //   overdampedMove(_network, m_dt);
  //   energy = _network.computeEnergy();
  // }
  // m_nextdt = std::clamp(2 * (m_prevEnergy - energy) / FdotF, dtMin, dtMax);
  // m_prevEnergy = energy;
}

auto networkV4::SteepestDescent::getDt() const -> double
{
  return m_dt;
}

auto networkV4::SteepestDescent::eGradient(network& _network,
                                           const double _step) const -> double
{
  ///const double norm = tools::norm(_network.getNodes().forces());
  overdampedMove(_network, - _step);
  const double energy = _network.computeEnergy();
  overdampedMove(_network, 2.0 * _step);
  const double newEnergy = _network.computeEnergy();
  overdampedMove(_network, -_step);
  return (newEnergy - energy) / ( 2.0 * _step);
}

auto networkV4::SteepestDescent::lineSearch(network& _network) const -> double
{
    //double a = config::intergrators::adaptiveIntergrator::dtMax;
    //double tau = 0.5;
    //double c = 1e-2;
//
    //double initialEnergy = _network.computeEnergy();
    //double fdotf = tools::xdoty(_network.getNodes().forces(), _network.getNodes().forces());
    //double Ea;
//
    //overdampedMove(_network, a);
    //Ea = _network.computeEnergy();
    //overdampedMove(_network, -a);
    //while (Ea > initialEnergy - a * c * fdotf) {
    //    a *= tau;
    //    overdampedMove(_network, a);
    //    Ea = _network.computeEnergy();
    //    overdampedMove(_network, -a);
    //    if (a < config::intergrators::adaptiveIntergrator::dtMin) {
    //        a = config::intergrators::adaptiveIntergrator::dtMin;
    //        break;
    //    }
    //}
    //return a;


  //m_h = m_nexth;
  //double energy;
//
  //std::size_t maxIter = config::intergrators::adaptiveIntergrator::maxIter;
  //double qMin = config::intergrators::adaptiveHeun::qMin;
  //double qMax = config::intergrators::adaptiveHeun::qMax;
  //double dtMin = config::intergrators::adaptiveIntergrator::dtMin;
  //double dtMax = config::intergrators::adaptiveIntergrator::dtMax;
//
  //std::size_t iters = 0;
//
  //nodes& networkNodes = _network.getNodes();
  //_network.computeForces();
  //if (std::isnan(m_prevEnergy)) {
  //  m_prevEnergy = _network.getEnergy();
  //}
  //const double maxLen = tools::maxLength(networkNodes.forces());
//
  //m_tempPositions = networkNodes.positions();
  //m_tempForces = networkNodes.forces();
//
  //while (iters < maxIter) {
  //  m_dt = std::clamp(m_h / maxLen, dtMin, dtMax);
  //  overdampedMove(_network, m_dt);
  //  energy = _network.computeEnergy();
//
  //  if (energy < m_prevEnergy) {
  //    m_prevEnergy = energy;
  //    m_nexth = m_h * 1.2;
  //    break;
  //  }
  //  networkNodes.positions() = m_tempPositions;
  //  networkNodes.forces() = m_tempForces;
  //  m_h *= 0.2;
  //  ++iters;
  //}
  //if (iters == maxIter) {
  //  std::runtime_error("Max iterations reached");
  //}

  double a = config::intergrators::adaptiveIntergrator::dtMin;
  double b = config::intergrators::adaptiveIntergrator::dtMax;
  double tol = 1e-8;

  double initialEnergy = _network.computeEnergy();
  // overdampedMove(_network, a);
  double fa = eGradient(_network);
  // double Ea = _network.computeEnergy();
  // overdampedMove(_network, -a);
  overdampedMove(_network, b);
  double fb = eGradient(_network);
  double Eb = _network.computeEnergy();
  overdampedMove(_network, -b);

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

    overdampedMove(_network, xITP);
    double fITP = eGradient(_network);
    overdampedMove(_network, -xITP);
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
    : m_dt(config::intergrators::default_dt)
    , m_tol(config::intergrators::miminizer::tol)
{
}

networkV4::FireMinimizer::FireMinimizer(double _tol)
    : m_dt(config::intergrators::default_dt)
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
  size_t maxIter = config::intergrators::miminizer::maxIter;
  double dtMin = config::intergrators::adaptiveIntergrator::dtMin;
  double dtMax = config::intergrators::adaptiveIntergrator::dtMax;

  double alpha0 = config::intergrators::fire2::alpha0;
  double finc = config::intergrators::fire2::finc;
  double fdec = config::intergrators::fire2::fdec;
  double falpha = config::intergrators::fire2::falpha;
  size_t Ndelay = config::intergrators::fire2::Ndelay;
  size_t Nnegmax = config::intergrators::fire2::Nnegmax;
  double dmax = config::intergrators::fire2::dmax;

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
  updateVelocities(_network, m_dt);

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
      move(_network, -0.5 * m_dt);
      networkNodes.clearVelocities();
      _network.computeForces();
    }
    updateVelocities(_network, m_dt);
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
    move(_network, m_dt);
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
  size_t maxIter = config::intergrators::miminizer::maxIter;

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
