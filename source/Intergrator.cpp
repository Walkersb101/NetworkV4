#include <algorithm>
#include <cstddef>
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

networkV4::intergrator::intergrator() {}

networkV4::intergrator::~intergrator() {}

void networkV4::intergrator::integrate(network& _network) {}

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
  double qMin = config::adaptiveIntergrator::qMin;
  double qMax = config::adaptiveIntergrator::qMax;
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