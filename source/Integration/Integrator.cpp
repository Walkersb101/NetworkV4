#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "Integrator.hpp"

#include "Core/Network.hpp"
#include "Core/Nodes.hpp"
#include "Core/Tensor2.hpp"
#include "Core/Vec2.hpp"
#include "Misc/Config.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"

networkV4::integrator::integrator() {}

networkV4::integrator::~integrator() {}

void networkV4::integrator::integrate(network& _network) {}

auto networkV4::integrator::getDt() const -> double
{
  return 0.0;
}

void networkV4::Integrationtools::overdampedMove(networkV4::network& _network,
                                                 double _dt)
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

void networkV4::Integrationtools::move(networkV4::network& _network, double _dt)
{
  networkV4::nodes& nodes = _network.getNodes();
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (nodes.fixed(i)) {
      continue;
    }
    nodes.position(i) += nodes.velocity(i) * _dt;
    _network.wrapPosition(nodes.position(i));
  }
}

void networkV4::Integrationtools::updateVelocities(networkV4::network& _network,
                                                   double _alpha)
{
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.velocity(i) += _nodes.force(i) * _alpha / _nodes.mass(i);
  }
}

void networkV4::Integrationtools::heunAverage(
    networkV4::network& _network,
    const std::vector<vec2d>& _tempForces,
    double _dt)
{
  networkV4::nodes& nodes = _network.getNodes();
  for (std::size_t i = 0; i < nodes.size(); ++i) {
    if (nodes.fixed(i)) {
      continue;
    }
    nodes.position(i) += (nodes.force(i) - _tempForces[i]) * _dt * 0.5;
    _network.wrapPosition(nodes.position(i));
    nodes.force(i) = (nodes.force(i) + _tempForces[i]) * 0.5;
  }
}
void networkV4::Integrationtools::copyVector(const std::vector<vec2d>& _src,
                                             std::vector<vec2d>& _dst)
{
  _dst.resize(_src.size());
  std::copy(_src.begin(), _src.end(), _dst.begin());
}

auto networkV4::Integrationtools::vectorDiffNorm(const std::vector<vec2d>& _src,
                                                 const std::vector<vec2d>& _dst,
                                                 const std::vector<bool>& _mask,
                                                 bool _invertmask)
    -> double
{
  if (_src.size() != _dst.size() || _src.size() != _mask.size()) {
    throw std::invalid_argument(
        "vectorDiffNorm: vectors must be of equal size");
  }
  double sum = 0.0;
  for (std::size_t i = 0; i < _src.size(); ++i) {
    sum += (_mask[i] != _invertmask) * (_src[i] - _dst[i]).lengthSquared();
  }
  return std::sqrt(sum);
}