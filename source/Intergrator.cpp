#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "Intergrator.hpp"

#include "Config.hpp"
#include "Network.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"

auto averageStresses(
    const std::unordered_map<networkV4::bondType, tensor2d>& _stresses1,
    const std::unordered_map<networkV4::bondType, tensor2d>& _stresses2)
    -> std::unordered_map<networkV4::bondType, tensor2d>
{
  // Check if _stresses2 has the same keys as _stresses1
  for (const auto& i : _stresses1) {
    if (_stresses2.find(i.first) == _stresses2.end()) {
      throw std::runtime_error(
          "stresses2 does not have the same keys as stresses1");
    }
  }
  std::unordered_map<networkV4::bondType, tensor2d> avgStresses = _stresses1;
  for (auto& i : avgStresses) {
    avgStresses[i.first] += _stresses2.at(i.first);
    avgStresses[i.first] *= 0.5;
  }
  return avgStresses;
}

void overdampedMove(networkV4::network& _network, double _dt)
{
  networkV4::nodes& _nodes = _network.getNodes();
#pragma omp parallel for
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
#pragma omp parallel for
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
#pragma omp parallel for
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
#pragma omp parallel for
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

auto norm(const std::vector<vec2d>& _vec) -> double
{
  double sum = std::accumulate(_vec.begin(),
                               _vec.end(),
                               0.0,
                               [](double _sum, const vec2d& _v)
                               { return _sum + _v.lengthSquared(); });
  return std::sqrt(sum);
}

auto maxAbsComponent(const std::vector<vec2d>& _vec) -> double
{
  return std::accumulate(_vec.begin(),
                         _vec.end(),
                         0.0,
                         [](double _max, const vec2d& _v)
                         { return std::max(_max, _v.abs().max()); });
  /*
    double max = 0.0;
    for (size_t i = 0; i < _vec.size(); ++i) {
      double val = _vec[i].abs().max();
      if (val > max) {
        max = val;
      }
    }
    return max;
  */
}

networkV4::intergrator::intergrator() {}

networkV4::intergrator::~intergrator() {}

void networkV4::intergrator::integrate(network& _network) {}

auto networkV4::intergrator::getDt() const -> double
{
  return 0.0;
}

networkV4::overdampedEuler::overdampedEuler()
    : m_dt(config::default_dt)
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
    : m_dt(config::default_dt)
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
      averageStresses(_network.getStresses(), m_tempStresses);
}

auto networkV4::OverdampedEulerHeun::getDt() const -> double
{
  return m_dt;
}

networkV4::OverdampedAdaptiveEulerHeun::OverdampedAdaptiveEulerHeun()
    : m_dt(config::default_dt)
    , m_nextdt(config::default_dt)
    , m_esp(config::adaptiveIntergrator::esp)
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

  std::size_t maxIter = config::adaptiveIntergrator::maxIter;
  double qMin = config::adaptiveHeun::qMin;
  double qMax = config::adaptiveHeun::qMax;
  double dtMin = config::adaptiveIntergrator::dtMin;
  double dtMax = config::adaptiveIntergrator::dtMax;

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

  heunAverage(_network, m_tempForces, m_dt);

  _network.getStresses() =
      averageStresses(_network.getStresses(), m_tempStresses);

  m_nextdt = std::clamp(m_dt * q, dtMin, dtMax);
}

auto networkV4::OverdampedAdaptiveEulerHeun::getDt() const -> double
{
  return m_dt;
}

auto networkV4::OverdampedAdaptiveEulerHeun::forceErrorNorm(nodes& _nodes)
    -> double
{
  double maxForce = 0.0;
#pragma omp parallel for reduction(std::max : max)
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    double force =
        _nodes.fixed(i) ? 0.0 : (m_tempForces[i] - _nodes.force(i)).abs().max();
    maxForce = force > maxForce ? force : maxForce;
  }
  return maxForce;
}

networkV4::FireMinimizer::FireMinimizer()
    : m_dt(config::default_dt)
    , m_tol(config::miminizer::tol)
{
}

networkV4::FireMinimizer::FireMinimizer(double _tol)
    : m_dt(config::default_dt)
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
  size_t maxIter = config::miminizer::maxIter;
  double dtMin = config::adaptiveIntergrator::dtMin;
  double dtMax = config::adaptiveIntergrator::dtMax;

  double alpha0 = config::fire2::alpha0;
  double finc = config::fire2::finc;
  double fdec = config::fire2::fdec;
  double falpha = config::fire2::falpha;
  size_t Ndelay = config::fire2::Ndelay;
  size_t Nnegmax = config::fire2::Nnegmax;

  size_t Npos = 0;
  size_t Nneg = 0;
  double alpha = alpha0;

  nodes& networkNodes = _network.getNodes();
  networkNodes.clearVelocities();
  _network.computeForces();

  if (norm(networkNodes.forces()) < m_tol) {
    return;
  }

  size_t iters = 0;
  while (iters < maxIter) {
    double P = power(_network);

    if (P >= 0.0) {
      Npos++;
      Nneg = 0;
      if (Npos > Ndelay) {
        m_dt = std::min(m_dt * finc, dtMax);
        alpha *= falpha;
      }
    } else {
      Npos = 0;
      Nneg++;
      if (Nneg > Nnegmax) {
        break;
      }
      if (iters > Ndelay) {
        m_dt = std::max(m_dt * fdec, dtMin);
        alpha = alpha0;
      }
      move(_network, -0.5 * m_dt);
      _network.computeForces();
      networkNodes.clearVelocities();
    }

    updateVelocities(_network, 0.5 * m_dt);
    double fScale =
        alpha * norm(networkNodes.velocities()) / norm(networkNodes.forces());
    double vScale = 1.0 - alpha;
    for (std::size_t i = 0; i < networkNodes.size(); ++i) {
      if (networkNodes.fixed(i)) {
        continue;
      }
      networkNodes.velocity(i) *= vScale;
      networkNodes.velocity(i) += networkNodes.force(i) * fScale;
    }
    move(_network, m_dt);
    _network.computeForces();
    updateVelocities(_network, 0.5 * m_dt);

    double error = maxAbsComponent(networkNodes.forces());
    if (error < m_tol) {
      break;
    }
    ++iters;
  }
}

auto networkV4::FireMinimizer::power(const network& _network) -> double
{
  const nodes& networkNodes = _network.getNodes();
  double power = 0.0;
  for (std::size_t i = 0; i < networkNodes.size(); ++i) {
    power += networkNodes.force(i).dot(networkNodes.velocity(i));
  }
  return power;
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
  size_t maxIter = config::miminizer::maxIter;

  size_t iters = 0;
  while (iters < maxIter) {
    m_integrator.integrate(_network);
    double error = maxAbsComponent(_network.getNodes().forces());
    if (error < m_tol) {
      break;
    }
    ++iters;
  }
}