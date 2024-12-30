#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Core/Network.hpp"
#include "Core/Nodes.hpp"

#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{


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
void overdampedMove(network& _network, double _dt, double _gamma = 1.0);
void move(network& _network, double _dt);
void updateVelocities(network& _network, double _alpha);
void heunAverageOverdamped(network& _network,
                 const std::vector<Utils::vec2d>& _tempForces,
                 double _dt,
                 double _gamma = 1.0);
void copyVector(const std::vector<Utils::vec2d>& _src, std::vector<Utils::vec2d>& _dst);
auto vectorDiffNorm(const std::vector<Utils::vec2d>& _src,
                    const std::vector<Utils::vec2d>& _dst) -> double;
}  // namespace Integrationtools

}  // namespace networkV4
