#include "Integrator.hpp"

#include <algorithm>
#include <cstddef>
#include <numeric>
#include <stdexcept>

#include <range/v3/all.hpp>

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

  double invgamma = 1.0 / _gamma;
  double overdampedScale = _dt * invgamma;

  std::transform(positions.begin(),
                 positions.end(),
                 _tempForces.begin(),
                 positions.begin(),
                 [overdampedScale](const auto& position, const auto& force)
                 { return position - force * overdampedScale; }); // x_{n} = x'_{n+1} - f(x_{n}) * dt / gamma
  std::transform(forces.begin(),
                 forces.end(),
                 _tempForces.begin(),
                 forces.begin(),
                 [](const auto& force, const auto& tempForce)
                 { return (force + tempForce) * 0.5; }); // f_{n+1} = 0.5 * (f(x_{n}) + f(x'_{n+1}))
  std::transform(forces.begin(),
                 forces.end(),
                 velocities.begin(),
                 [invgamma](const auto& force)
                 { return force * invgamma; }); // v_{n+1} = f_{n+1} / gamma
  std::transform(positions.begin(),
                 positions.end(),
                 velocities.begin(),
                 positions.begin(),
                 [_dt](const auto& position, const auto& velocity)
                 { return position + velocity * _dt; }); // x_{n+1} = x_{n} + v_{n+1} * dt
}

void networkV4::Integrationtools::copyVector(
    const std::vector<Utils::vec2d>& _src, std::vector<Utils::vec2d>& _dst)
{
  _dst.resize(_src.size());
  std::copy(_src.begin(), _src.end(), _dst.begin());
}

#include <range/v3/all.hpp>

auto networkV4::Integrationtools::vectorDiffNorm(const std::vector<Utils::vec2d>& _src,
                                                 const std::vector<Utils::vec2d>& _dst,
                                                 const std::vector<bool>& _mask,
                                                 bool _invertmask) -> double
{
  if (_src.size() != _dst.size() || _src.size() != _mask.size()) {
    throw std::invalid_argument(
        "vectorDiffNorm: vectors must be of equal size");
  }

  auto sum = ranges::accumulate(
      ranges::views::zip(_src, _dst, _mask) |
      ranges::views::transform([_invertmask](const auto& t) {
        const auto& [src, dst, mask] = t;
        return mask != _invertmask ? (src - dst).normSquared() : 0.0;
      }),
      0.0);

  return std::sqrt(sum);
}