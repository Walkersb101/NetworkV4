#pragma once

#include <unordered_map>

#include "Bonds.hpp"
#include "Intergrator.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"

namespace network
{
class network
{
public:
  network();
  ~network();

public:
    auto getNodes() -> nodes&;
    auto getNodes() const -> const nodes&;

    auto getStresses() -> std::unordered_map<bondType, tensor2d>&;

public:
  void computeForces();

private:
  auto minDist(const vec2d& _pos1, const vec2d& _pos2) const -> vec2d;

public:
    void clearStresses();
    void normalizeStresses();

private:
  vec2d m_restSize;
  vec2d m_domain;

  nodes m_nodes;
  bonds m_bonds;

  std::unordered_map<bondType, tensor2d> m_stresses;

  double m_shearStrain;
  double m_shearStraingRate;

  double m_elongationStrain;
  double m_elongationStraingRate;
};
}  // namespace network