#include <algorithm>

#include "Nodes.hpp"

networkV4::nodes::nodes()
    : m_count(0)
{
}

void networkV4::nodes::clear()
{
  m_count = 0;
  m_positions.clear();
  m_velocities.clear();
  m_forces.clear();
  m_masses.clear();
}

void networkV4::nodes::resize(std::size_t _count)
{
  m_count = _count;
  m_positions.resize(_count);
  m_velocities.resize(_count);
  m_forces.resize(_count);
  m_masses.resize(_count);
  m_fixed.resize(_count);
}

void networkV4::nodes::reserve(std::size_t _count)
{
  m_positions.reserve(_count);
  m_velocities.reserve(_count);
  m_forces.reserve(_count);
  m_masses.reserve(_count);
  m_fixed.reserve(_count);
}

void networkV4::nodes::shrink_to_fit()
{
  m_positions.shrink_to_fit();
  m_velocities.shrink_to_fit();
  m_forces.shrink_to_fit();
  m_masses.shrink_to_fit();
  m_fixed.shrink_to_fit();
}

auto networkV4::nodes::size() const -> std::size_t
{
  return m_count;
}

auto networkV4::nodes::empty() const -> bool
{
  return m_count == 0;
}

auto networkV4::nodes::position(std::size_t _idx) const -> vec2d
{
  return m_positions[_idx];
}

auto networkV4::nodes::position(std::size_t _idx) -> vec2d&
{
  return m_positions[_idx];
}

auto networkV4::nodes::velocity(std::size_t _idx) const -> vec2d
{
  return m_velocities[_idx];
}

auto networkV4::nodes::velocity(std::size_t _idx) -> vec2d&
{
  return m_velocities[_idx];
}

auto networkV4::nodes::force(std::size_t _idx) const -> vec2d
{
  return m_forces[_idx];
}

auto networkV4::nodes::force(std::size_t _idx) -> vec2d&
{
  return m_forces[_idx];
}

auto networkV4::nodes::mass(std::size_t _idx) const -> double
{
  return m_masses[_idx];
}

auto networkV4::nodes::mass(std::size_t _idx) -> double&
{
  return m_masses[_idx];
}

auto networkV4::nodes::fixed(std::size_t _idx) const -> bool
{
  return m_fixed[_idx];
}

void networkV4::nodes::fix(std::size_t _idx)
{
  m_fixed[_idx] = true;
}

void networkV4::nodes::unfix(std::size_t _idx)
{
  m_fixed[_idx] = false;
}

auto networkV4::nodes::positions() -> std::vector<vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() -> std::vector<vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() -> std::vector<vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() -> std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::fixed() -> std::vector<bool>&
{
  return m_fixed;
}

auto networkV4::nodes::positions() const -> const std::vector<vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() const -> const std::vector<vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() const -> const std::vector<vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() const -> const std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::fixed() const -> const std::vector<bool>&
{
  return m_fixed;
}

void networkV4::nodes::addNode(const vec2d& _position)
{
  addNode(_position, vec2d(0.0, 0.0), vec2d(0.0, 0.0), 1.0);
}

void networkV4::nodes::addNode(const vec2d& _position,
                               const vec2d& _velocity,
                               const vec2d& _force,
                               double _mass)
{
  m_positions.push_back(_position);
  m_velocities.push_back(_velocity);
  m_forces.push_back(_force);
  m_masses.push_back(_mass);
  m_fixed.push_back(false);
  m_count++;
}

void networkV4::nodes::addNode(const vec2d& _position,
                               const vec2d& _velocity,
                               const vec2d& _force,
                               double _mass,
                               bool _fixed)
{
  addNode(_position, _velocity, _force, _mass);
  m_fixed.back() = _fixed;
}

void networkV4::nodes::removeNode(std::size_t _idx)
{
  m_positions.erase(m_positions.begin() + _idx);
  m_velocities.erase(m_velocities.begin() + _idx);
  m_forces.erase(m_forces.begin() + _idx);
  m_masses.erase(m_masses.begin() + _idx);
  m_fixed.erase(m_fixed.begin() + _idx);
  m_count--;
}

void networkV4::nodes::setNode(std::size_t _idx,
                               const vec2d& _position,
                               const vec2d& _velocity,
                               const vec2d& _force,
                               double _mass)
{
  m_positions[_idx] = _position;
  m_velocities[_idx] = _velocity;
  m_forces[_idx] = _force;
  m_masses[_idx] = _mass;
}

void networkV4::nodes::setNode(std::size_t _idx,
                               const vec2d& _position,
                               const vec2d& _velocity,
                               const vec2d& _force,
                               double _mass,
                               bool _fixed)
{
  setNode(_idx, _position, _velocity, _force, _mass);
  m_fixed[_idx] = _fixed;
}

void networkV4::nodes::clearVelocities()
{
  std::fill(PAR m_velocities.begin(), m_velocities.end(), vec2d(0.0, 0.0));
}

void networkV4::nodes::clearForces()
{
  std::fill(PAR m_forces.begin(), m_forces.end(), vec2d(0.0, 0.0));
}

void networkV4::nodes::zeroFixedForces()
{
  std::transform(PAR m_fixed.begin(),
                 m_fixed.end(),
                 m_forces.begin(),
                 m_forces.begin(),
                 [](bool _fixed, const vec2d& _force)
                 { return _fixed ? vec2d(0.0, 0.0) : _force; });
}