#pragma once

#include <algorithm>
#include <vector>

#include <range/v3/view/zip.hpp>

#include "Adaptive.hpp"
#include "Core/Network.hpp"
#include "Overdamped.hpp"

namespace networkV4
{
namespace integration
{

class AdaptiveOverdampedEulerHeun
    : public Overdamped
    , public Adaptive
{
public:
  AdaptiveOverdampedEulerHeun() = default;
  AdaptiveOverdampedEulerHeun(double _gamma, const AdaptiveParams& _params)
      : Overdamped(_gamma)
      , Adaptive(_params)
  {
  }
  ~AdaptiveOverdampedEulerHeun() override = default;

  auto step(network& _network, double _dt) -> dtOut
  {
    const double overdampedScale = _dt * m_invGamma;

    networkV4::nodes& nodes = _network.getNodes();
    m_fn = nodes.forces();
    std::transform(nodes.positions().begin(),
                   nodes.positions().end(),
                   nodes.forces().begin(),
                   nodes.positions().begin(),
                   [_dt](const auto& position, const auto& force)
                   { return position + force * overdampedScale; });

    _network.computeForces();

    const double halfOverdampedScale = 0.5 * overdampedScale;
    for (auto [position, velocity, force, tempForce] : ranges::views::zip(
             nodes.positions(), nodes.velocities(), nodes.forces(), m_fn))
    {
      position += (force - tempForce) * halfOverdampedScale;
      force = (force + tempForce) * 0.5;
      velocity = force * m_invGamma;
    }
  }

private:
  std::vector<Utils::vec2d> m_fn;
};

}  // namespace integration
}  // namespace networkV4
