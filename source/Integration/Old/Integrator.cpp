#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include "Integrator.hpp"

#include <range/v3/numeric/accumulate.hpp>
#include <range/v3/range/conversion.hpp>
#include <range/v3/view/transform.hpp>
#include <range/v3/view/zip.hpp>

#include "Core/Network.hpp"
#include "Core/Nodes.hpp"
#include "Misc/Config.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Tools.hpp"
#include "Misc/Vec2.hpp"

networkV4::integrator::integrator() {}

networkV4::integrator::~integrator() {}

void networkV4::integrator::integrate(network& _network) {}

auto networkV4::integrator::getDt() const -> double
{
  return 0.0;
}

void networkV4::Integrationtools::overdampedMove(networkV4::network& _network,
                                                 double _dt,
                                                 double _gamma)
{
  auto& positions = _network.getNodes().positions();
  auto& velocities = _network.getNodes().velocities();
  const auto& forces = _network.getNodes().forces();

  double invGamma = 1.0 / _gamma;

  std::transform(forces.begin(),
                 forces.end(),
                 velocities.begin(),
                 [invGamma](const auto& force) { return force * invGamma; });
  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 positions.begin(),
                 [_dt](const auto& position, const auto& velocity)
                 { return position + velocity * _dt; });
}

void networkV4::Integrationtools::move(networkV4::network& _network, double _dt)
{
  auto& positions = _network.getNodes().positions();
  const auto& velocities = _network.getNodes().velocities();

  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 positions.begin(),
                 [_dt](const auto& position, const auto& velocity)
                 { return position + velocity * _dt; });
}

void networkV4::Integrationtools::updateVelocities(networkV4::network& _network,
                                                   double _alpha)
{
  auto& velocities = _network.getNodes().velocities();
  const auto& forces = _network.getNodes().forces();
  std::transform(velocities.begin(),
                 velocities.end(),
                 forces.begin(),
                 velocities.begin(),
                 [_alpha](const auto& velocity, const auto& force)
                 { return velocity + force * _alpha; });
}

void networkV4::Integrationtools::heunAverageOverdamped(
    networkV4::network& _network,
    const std::vector<Utils::vec2d>& _tempForces,
    double _dt,
    double _gamma)
{
  auto& positions = _network.getNodes().positions();
  auto& velocities = _network.getNodes().velocities();
  auto& forces = _network.getNodes().forces();

  double overdampedScale = _dt / _gamma;

  for (auto [position, velocity, force, tempForce] :
       ranges::views::zip(positions, velocities, forces, _tempForces))
  {
    position += (force - tempForce) * overdampedScale * 0.5;
    force = (force + tempForce) * 0.5;
  }
}

void networkV4::Integrationtools::copyVector(
    const std::vector<Utils::vec2d>& _src, std::vector<Utils::vec2d>& _dst)
{
  _dst.resize(_src.size());
  std::copy(_src.begin(), _src.end(), _dst.begin());
}

auto networkV4::Integrationtools::vectorDiffNorm(
    const std::vector<Utils::vec2d>& _src,
    const std::vector<Utils::vec2d>& _dst) -> double
{
  if (_src.size() != _dst.size()) {
    throw std::invalid_argument(
        "vectorDiffNorm: vectors must be of equal size");
  }
  double sum = 0;
  for (auto [src, dst] : ranges::views::zip(_src, _dst)) {
    sum += (src - dst).normSquared();
  }
  return std::sqrt(sum);
}