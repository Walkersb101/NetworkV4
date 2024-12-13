#pragma once

#include <filesystem>
#include <vector>

#include "Misc/EnumMap.hpp"
#include "Misc/Enums.hpp"
#include "Core/Nodes.hpp"
#include "Core/Bonds.hpp"
#include "Core/box.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{

using stressMap =
    EnumMap<bondType, tensor2d, bondType::single, bondType::matrix>;

class network
{
public:

public:
  auto getNodes() -> nodes&;
  auto getNodes() const -> const nodes&;

  auto getBonds() -> bonds&;
  auto getBonds() const -> const bonds&;

  auto getStresses() -> stressMap&;
  auto getStresses() const -> const stressMap&;

  auto getEnergy() const -> const double;

  auto getShearStrain() const -> const double;
  auto getElongationStrain() const -> const double;

  auto getRestSize() const -> const vec2d;
  auto getDomain() const -> const vec2d;

  void setDomain(const vec2d& _domain);

public:
  void computeForces();

  auto computeEnergy() -> double;

  auto minDist(const vec2d& _pos1, const vec2d& _pos2) const -> vec2d;
  void wrapPosition(vec2d& _pos) const;
  auto wrapedPosition(const vec2d& _pos) const -> vec2d;

public:
  void shear(double _step);
  void elongate(double _step);

private:
  void applyBond(const bond& _bond,
                 const nodes& _nodes,
                 std::vector<vec2d>& _forces,
                 stressMap& _stresses,
                 double& _energy);
  void applyBond(const bond& _bond);

public:
  auto bondStrain(const bond& _bond) const -> double;
  auto bondEnergy(const bond& _bond) const -> double;

public:
  void initStresses();
  void clearStresses();
  void normalizeStresses();
  auto getGlobalStress() const -> tensor2d;

private:
  box m_domain;
  box m_restSize;

  nodes m_nodes;
  bonded::bonds m_bonds;

  double m_energy;
};

}  // namespace networkV4
