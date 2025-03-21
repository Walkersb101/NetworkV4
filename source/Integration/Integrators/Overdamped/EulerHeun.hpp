#pragma once

#include <algorithm>
#include <vector>

#include <range/v3/view/zip.hpp>

#include "Core/Network.hpp"
#include "Integration/Integrators/IntegratorBase.hpp"
#include "OverdampedBase.hpp"

namespace networkV4
{
namespace integration
{

class OverdampedEulerHeun
    : public Overdamped
    , public IntegratorBase
{
public:
  OverdampedEulerHeun() = default;
  OverdampedEulerHeun(double _zeta, double _dt)
      : Overdamped(_zeta)
      , IntegratorBase(_dt)
  {
  }
  ~OverdampedEulerHeun() override = default;

  void step(network& _network) override
  {
    const double overdampedScale = m_dt * m_invZeta;

    networkV4::nodes& nodes = _network.getNodes();
    m_fn = nodes.forces();
    std::transform(nodes.positions().begin(),
                   nodes.positions().end(),
                   nodes.forces().begin(),
                   nodes.positions().begin(),
                   [overdampedScale](const auto& position, const auto& force)
                   { return position + force * overdampedScale; });

    _network.computeForces();

    const double halfOverdampedScale = 0.5 * overdampedScale;
    for (auto [position, velocity, force, tempForce] : ranges::views::zip(
             nodes.positions(), nodes.velocities(), nodes.forces(), m_fn))
    {
      position += (force - tempForce) * halfOverdampedScale;
      force = (force + tempForce) * 0.5;
      velocity = force * m_invZeta;
    }
  }

private:
  std::vector<Utils::Math::vec2d> m_fn;
};

}  // namespace integration
}  // namespace networkV4
