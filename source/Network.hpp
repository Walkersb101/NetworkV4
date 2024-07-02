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

private:
  void loadFromBinV1(const std::filesystem::path& _filename, double _lambda);
  void loadFromBinV2(const std::filesystem::path& _filename);

public:
  void loadFromBin(const std::filesystem::path& _filename,
                   loadVersion _version,
                   double _lambda = 1.0);

public:
  auto getNodes() -> nodes&;
  auto getNodes() const -> const nodes&;

  auto getStresses() -> std::unordered_map<bondType, tensor2d>&;

public:
  void computeForces();

  auto minDist(const vec2d& _pos1, const vec2d& _pos2) const -> vec2d;
  void wrapPosition(vec2d& _pos) const;
  auto wrapedPosition(const vec2d& _pos) const -> vec2d;

private:
  void applyBond(const bond& _bond,
                 const nodes& _nodes,
                 std::vector<vec2d>& _forces,
                 std::unordered_map<bondType, tensor2d>& _stresses);
  void applyBond(const bond& _bond);

public:
  void clearStresses();
  void normalizeStresses();

private:
  vec2d m_restSize;
  vec2d m_domain;

  nodes m_nodes;
  bonds m_bonds;

  std::unordered_map<bondType, tensor2d> m_stresses;
  tensor2d m_avgStress;

  double m_shearStrain;
  double m_shearStraingRate;

  double m_elongationStrain;
  double m_elongationStraingRate;
};
}  // namespace networkV4

//#pragma omp declare reduction( \
//        stresses_plus : std::unordered_map<bondType, tensor2d> : std:: \
//        transform(omp_out.begin(), \
//                      omp_out.end(), \
//                      omp_in.begin(), \
//                      omp_out.begin(), \
//                      std::plus<tensor2d>())) initializer(omp_priv = omp_orig)