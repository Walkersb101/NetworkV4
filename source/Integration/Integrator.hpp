#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Core/Network.hpp"
#include "Core/Nodes.hpp"
#include "Core/Tensor2.hpp"
#include "Core/Vec2.hpp"
#include "Misc/EnumMap.hpp"

namespace networkV4
{

using stressMap = EnumMap<bondType, tensor2d, bondType::single, bondType::matrix>;

class integrator
{
public:
  integrator();

  virtual ~integrator();

public:
  virtual void integrate(network& _network);
  virtual auto getDt() const -> double;
};

namespace Integrationtools
{
    void overdampedMove(network& _network, double _dt);
    void move(network& _network, double _dt);
    void updateVelocities(network& _network, double _alpha);
    void heunAverage(network& _network, const std::vector<vec2d>& _tempForces, double _dt);
}  // namespace itergrationTools

}  // namespace networkV4
