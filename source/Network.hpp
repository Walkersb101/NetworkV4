#pragma once

#include <filesystem>
#include <unordered_map>
#include <vector>

#include "Bonds.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"

namespace networkV4
{
enum loadVersion : std::uint8_t
{
  binV1 = 1,
  binV2 = 2
};

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

  auto getStresses() -> std::unordered_map<bondType, tensor2d>&;
  auto getStresses() const -> const std::unordered_map<bondType, tensor2d>&;

  auto getShearStrain() const -> double;
  auto getElongationStrain() const -> double;

  auto getShearStrain() -> double&;
  auto getElongationStrain() -> double&;

public:
  void computeForces();

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
                 std::unordered_map<bondType, tensor2d>& _stresses);
  void applyBond(const bond& _bond);

public:
  auto bondStrain(const bond& _bond) const -> double;

public:
  void strainBreakData(double& _maxDistAbove, std::size_t& _count) const;
  auto strainBreak() -> std::vector<std::size_t>;

public:
  void initStresses();
  void clearStresses();
  void normalizeStresses();
  auto globalStress() const -> tensor2d;

private:
  vec2d m_restSize;
  vec2d m_domain;

  nodes m_nodes;
  bonds m_bonds;

  std::unordered_map<bondType, tensor2d> m_stresses;
  tensor2d m_avgStress;

  double m_shearStrain;

  double m_elongationStrain;
};

}  // namespace networkV4

//#pragma omp declare reduction( \
//        stresses_plus : std::unordered_map<bondType, tensor2d> : std:: \
//        transform(omp_out.begin(), \
//                      omp_out.end(), \
//                      omp_in.begin(), \
//                      omp_out.begin(), \
//                      std::plus<tensor2d>())) initializer(omp_priv = omp_orig)