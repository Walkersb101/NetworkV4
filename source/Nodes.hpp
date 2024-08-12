#pragma once

#include <cstdint>
#include <ranges>
#include <vector>

#include "Config.hpp"
#include "Vec2.hpp"

#pragma omp declare reduction( \
        vec_vec2d_plus : std::vector<vec2d> : std::transform(omp_out.begin(), \
                    omp_out.end(), \
                    omp_in.begin(), \
                    omp_out.begin(),  \
                    std::plus<vec2d>())) \
    initializer(omp_priv = decltype(omp_orig)(omp_orig.size()))

namespace networkV4
{
class nodes
{
public:
  nodes();

public:
  void clear();
  void resize(std::size_t _count);
  void reserve(std::size_t _count);
  void shrink_to_fit();

  auto size() const -> std::size_t;
  auto empty() const -> bool;

public:
  auto position(std::size_t _idx) const -> vec2d;
  auto position(std::size_t _idx) -> vec2d&;

  auto velocity(std::size_t _idx) const -> vec2d;
  auto velocity(std::size_t _idx) -> vec2d&;

  auto force(std::size_t _idx) const -> vec2d;
  auto force(std::size_t _idx) -> vec2d&;

  auto mass(std::size_t _idx) const -> double;
  auto mass(std::size_t _idx) -> double&;

  auto fixed(std::size_t _idx) const -> bool;
  void fix(std::size_t _idx);
  void unfix(std::size_t _idx);

public:
  auto positions() -> std::vector<vec2d>&;
  auto velocities() -> std::vector<vec2d>&;
  auto forces() -> std::vector<vec2d>&;
  auto masses() -> std::vector<double>&;
  auto fixed() -> std::vector<bool>&;

  auto positions() const -> const std::vector<vec2d>&;
  auto velocities() const -> const std::vector<vec2d>&;
  auto forces() const -> const std::vector<vec2d>&;
  auto masses() const -> const std::vector<double>&;
  auto fixed() const -> const std::vector<bool>&;

public:
  void addNode(const vec2d& _position);
  void addNode(const vec2d& _position,
               const vec2d& _velocity,
               const vec2d& _force,
               double _mass);
  void addNode(const vec2d& _position,
               const vec2d& _velocity,
               const vec2d& _force,
               double _mass,
               bool _fixed);
  void removeNode(std::size_t _idx);
  void setNode(std::size_t _idx,
               const vec2d& _position,
               const vec2d& _velocity,
               const vec2d& _force,
               double _mass);
  void setNode(std::size_t _idx,
               const vec2d& _position,
               const vec2d& _velocity,
               const vec2d& _force,
               double _mass,
               bool _fixed);

public:
  void clearVelocities();
  void clearForces();
  void zeroFixedForces();

private:
  std::size_t m_count;
  std::vector<vec2d> m_positions;
  std::vector<vec2d> m_velocities;
  std::vector<vec2d> m_forces;
  std::vector<double> m_masses;
  std::vector<bool> m_fixed;
};
}  // namespace networkV4