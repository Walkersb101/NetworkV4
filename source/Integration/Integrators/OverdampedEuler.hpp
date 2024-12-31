#pragma once

#include <algorithm>
#include <vector>

#include "Core/Network.hpp"
#include "IntegratorBase.hpp"
#include "Overdamped.hpp"

namespace networkV4
{
namespace integration
{

class OverdampedEuler
    : public Overdamped
    , public IntegratorBase
{
public:
  OverdampedEuler() = default;
  OverdampedEuler(double _gamma, double _dt)
      : Overdamped(_gamma)
      , IntegratorBase(_dt)
  {
  }
  ~OverdampedEuler() override = default;

  void step(network& _network) override
  {
    auto& nodes = _network.getNodes();

    std::transform(nodes.forces().begin(),
                   nodes.forces().end(),
                   nodes.velocities().begin(),
                   [this](const auto& force) { return force * m_invGamma; });
    std::transform(nodes.positions().begin(),
                   nodes.positions().end(),
                   nodes.velocities().begin(),
                   nodes.positions().begin(),
                   [this](const auto& position, const auto& velocity)
                   { return position + velocity * m_dt; });
  }
};

}  // namespace integration
}  // namespace networkV4
