#pragma once

#include <cstdint>
#include <vector>

#include "Config.hpp"
#include "Vec2.hpp"

namespace network
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
    auto positions_begin() -> std::vector<vec2d>::iterator;
    auto positions_end() -> std::vector<vec2d>::iterator;

    auto velocities_begin() -> std::vector<vec2d>::iterator;
    auto velocities_end() -> std::vector<vec2d>::iterator;

    auto forces_begin() -> std::vector<vec2d>::iterator;
    auto forces_end() -> std::vector<vec2d>::iterator;

    auto masses_begin() -> std::vector<double>::iterator;
    auto masses_end() -> std::vector<double>::iterator;

    auto fixed_begin() -> std::vector<bool>::iterator;
    auto fixed_end() -> std::vector<bool>::iterator;

public:
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

private:
  std::size_t m_count;
  std::vector<vec2d> m_positions;
  std::vector<vec2d> m_velocities;
  std::vector<vec2d> m_forces;
  std::vector<double> m_masses;
  std::vector<bool> m_fixed;
};
}  // namespace network