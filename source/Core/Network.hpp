#pragma once

#include <filesystem>
#include <vector>

#include "Misc/EnumMap.hpp"
#include "Misc/Enums.hpp"
#include "Core/Nodes.hpp"
#include "Core/Bonds.hpp"
#include "Core/Tensor2.hpp"
#include "Core/Vec2.hpp"

namespace networkV4
{

using stressMap =
    EnumMap<bondType, tensor2d, bondType::single, bondType::matrix>;

class network
{
public:
  network();
  ~network();

public:
  void loadFromBinV1(std::ifstream& _file, double _lambda);
  void loadFromBinV2(std::ifstream& _file);

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
  vec2d m_domain;
  vec2d m_restSize;

  vec2d m_halfDomain;
  vec2d m_yshift;

  nodes m_nodes;
  bonds m_bonds;

  stressMap m_stresses;
  tensor2d m_avgStress;

  double m_energy;

  double m_shearStrain;
  double m_elongationStrain;
};

}  // namespace networkV4
