#pragma once

#include <algorithm>
#include <vector>

#include <range/v3/view/zip.hpp>

#include "Core/Network.hpp"
#include "IntegratorBase.hpp"
#include "Overdamped.hpp"

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
  OverdampedEulerHeun(double _gamma, double _dt)
      : Overdamped(_gamma)
      , IntegratorBase(_dt)
  {
  }
  ~OverdampedEulerHeun() override = default;

  void step(network& _network) override
  {
    const double overdampedScale = m_dt * m_invGamma;

    networkV4::nodes& nodes = _network.getNodes();
    m_fn = nodes.forces();
    std::transform(nodes.positions().begin(),
                   nodes.positions().end(),
                   nodes.forces().begin(),
                   nodes.positions().begin(),
                   [](const auto& position, const auto& force)
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
