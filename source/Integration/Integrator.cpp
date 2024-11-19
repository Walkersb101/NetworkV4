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
#include "Misc/Tools.hpp"
#include "Misc/Roots.hpp"

void networkV4::Integrationtools::overdampedMove(networkV4::network& _network, double _dt)
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
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.position(i) += _nodes.velocity(i) * _dt;
    _network.wrapPosition(_nodes.position(i));
  }
}

void networkV4::Integrationtools::updateVelocities(networkV4::network& _network, double _alpha)
{
  networkV4::nodes& _nodes = _network.getNodes();
  for (std::size_t i = 0; i < _nodes.size(); ++i) {
    if (_nodes.fixed(i)) {
      continue;
    }
    _nodes.velocity(i) += _nodes.force(i) * _alpha / _nodes.mass(i);
  }
}

void networkV4::Integrationtools::heunAverage(networkV4::network& _network,
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

networkV4::integrator::integrator() {}

networkV4::integrator::~integrator() {}

void networkV4::integrator::integrate(network& _network) {}

auto networkV4::integrator::getDt() const -> double
{
  return 0.0;
}

