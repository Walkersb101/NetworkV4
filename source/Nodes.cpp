#include "Nodes.hpp"

network::nodes::nodes()
    : m_count(0)
{
}

void network::nodes::clear()
{
  m_count = 0;
  m_positions.clear();
  m_velocities.clear();
  m_forces.clear();
  m_masses.clear();
}

void network::nodes::resize(std::size_t _count)
{
  m_count = _count;
  m_positions.resize(_count);
  m_velocities.resize(_count);
  m_forces.resize(_count);
  m_masses.resize(_count);
}

void network::nodes::reserve(std::size_t _count)
{
  m_positions.reserve(_count);
  m_velocities.reserve(_count);
  m_forces.reserve(_count);
  m_masses.reserve(_count);
}

void network::nodes::shrink_to_fit()
{
  m_positions.shrink_to_fit();
  m_velocities.shrink_to_fit();
  m_forces.shrink_to_fit();
  m_masses.shrink_to_fit();
}

auto network::nodes::size() const -> std::size_t
{
  return m_count;
}

auto network::nodes::empty() const -> bool
{
  return m_count == 0;
}

auto network::nodes::position(std::size_t _idx) const -> vec2d
{
  return m_positions[_idx];
}

auto network::nodes::position(std::size_t _idx) -> vec2d&
{
  return m_positions[_idx];
}

auto network::nodes::velocity(std::size_t _idx) const -> vec2d
{
  return m_velocities[_idx];
}

auto network::nodes::velocity(std::size_t _idx) -> vec2d&
{
  return m_velocities[_idx];
}

auto network::nodes::force(std::size_t _idx) const -> vec2d
{
  return m_forces[_idx];
}

auto network::nodes::force(std::size_t _idx) -> vec2d&
{
  return m_forces[_idx];
}

auto network::nodes::mass(std::size_t _idx) const -> double
{
  return m_masses[_idx];
}

auto network::nodes::mass(std::size_t _idx) -> double&
{
  return m_masses[_idx];
}

auto network::nodes::fixed(std::size_t _idx) const -> bool
{
  return m_fixed[_idx];
}

void network::nodes::fix(std::size_t _idx)
{
  m_fixed[_idx] = true;
}

void network::nodes::unfix(std::size_t _idx)
{
  m_fixed[_idx] = false;
}

auto network::nodes::positions_begin() -> std::vector<vec2d>::iterator
{
  return m_positions.begin();
}

auto network::nodes::positions_end() -> std::vector<vec2d>::iterator
{
  return m_positions.end();
}

auto network::nodes::velocities_begin() -> std::vector<vec2d>::iterator
{
  return m_velocities.begin();
}

auto network::nodes::velocities_end() -> std::vector<vec2d>::iterator
{
  return m_velocities.end();
}

auto network::nodes::forces_begin() -> std::vector<vec2d>::iterator
{
  return m_forces.begin();
}

auto network::nodes::forces_end() -> std::vector<vec2d>::iterator
{
  return m_forces.end();
}

auto network::nodes::masses_begin() -> std::vector<double>::iterator
{
  return m_masses.begin();
}

auto network::nodes::masses_end() -> std::vector<double>::iterator
{
  return m_masses.end();
}

auto network::nodes::fixed_begin() -> std::vector<bool>::iterator
{
  return m_fixed.begin();
}

auto network::nodes::fixed_end() -> std::vector<bool>::iterator
{
  return m_fixed.end();
}

void network::nodes::addNode(const vec2d& _position,
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

void network::nodes::addNode(const vec2d& _position,
                              const vec2d& _velocity,
                              const vec2d& _force,
                              double _mass,
                              bool _fixed)
{
  addNode(_position, _velocity, _force, _mass);
  m_fixed.back() = _fixed;
}

void network::nodes::removeNode(std::size_t _idx)
{
  m_positions.erase(m_positions.begin() + _idx);
  m_velocities.erase(m_velocities.begin() + _idx);
  m_forces.erase(m_forces.begin() + _idx);
  m_masses.erase(m_masses.begin() + _idx);
  m_fixed.erase(m_fixed.begin() + _idx);
  m_count--;
}

void network::nodes::setNode(std::size_t _idx,
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

void network::nodes::setNode(std::size_t _idx,
             const vec2d& _position,
             const vec2d& _velocity,
             const vec2d& _force,
             double _mass,
             bool _fixed)
{
  setNode(_idx, _position, _velocity, _force, _mass);
  m_fixed[_idx] = _fixed;
}

void network::nodes::clearVelocities()
{
  std::fill(PAR m_velocities.begin(), m_velocities.end(), vec2d(0.0, 0.0));
}

void network::nodes::clearForces()
{
  std::fill(PAR m_forces.begin(), m_forces.end(), vec2d(0.0, 0.0));
}
