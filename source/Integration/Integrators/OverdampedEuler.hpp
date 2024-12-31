#pragma once

#include <algorithm>
#include <vector>

#include "Core/Nodes.hpp"
#include "Overdamped.hpp"

namespace networkV4
{
namespace integration
{

class OverdampedEuler : public Overdamped
{
public:
  OverdampedEuler() = default;
  OverdampedEuler(double _gamma)
      : Overdamped(_gamma)
  {
  }
  ~OverdampedEuler() override = default;

  void step(networkV4::nodes& _nodes, double _dt)
  {
    std::transform(_nodes.forces().begin(),
                   _nodes.forces().end(),
                   _nodes.velocities().begin(),
                   [this](const auto& force) { return force * m_invGamma; });
    std::transform(_nodes.positions().begin(),
                   _nodes.positions().end(),
                   _nodes.velocities().begin(),
                   _nodes.positions().begin(),
                   [_dt](const auto& position, const auto& velocity)
                   { return position + velocity * _dt; });
  }
};

}  // namespace integration
}  // namespace networkV4
