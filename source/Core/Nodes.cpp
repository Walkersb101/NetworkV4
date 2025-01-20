#include <algorithm>
#include <stdexcept>

#include "Core/Nodes.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>

#include "Misc/Utils.hpp"

networkV4::nodes::nodes(const size_t _size)
{
  reserve(_size);
}

void networkV4::nodes::clear()
{
  m_globalIndices.clear();

  m_positions.clear();
  m_forces.clear();
  m_velocities.clear();

  m_masses.clear();
}

void networkV4::nodes::reserve(size_t _size)
{
  m_globalIndices.reserve(_size);

  m_positions.reserve(_size);
  m_forces.reserve(_size);
  m_velocities.reserve(_size);

  m_masses.reserve(_size);
}

auto networkV4::nodes::size() const -> size_t
{
  if (m_globalIndices.size() != m_positions.size()
      || m_globalIndices.size() != m_velocities.size()
      || m_globalIndices.size() != m_masses.size()
      || m_globalIndices.size() != m_forces.size())
  {
    throw std::runtime_error("nodes::size: inconsistent sizes");
  }

  return m_globalIndices.size();
}

auto networkV4::nodes::addNode(const Utils::Math::vec2d& _position,
                               const Utils::Math::vec2d& _velocity,
                               double _mass) -> size_t
{
  size_t index = size();
  // TODO: check if tags have been set
  pushNode(index, _position, _velocity, _mass, Utils::Math::vec2d({0.0, 0.0}));
  // m_nextIndex++;
  return index;
}

auto networkV4::nodes::indices() const -> const std::vector<size_t>&
{
  return m_globalIndices;
}

auto networkV4::nodes::positions() const -> const std::vector<Utils::Math::vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() const -> const std::vector<Utils::Math::vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() const -> const std::vector<Utils::Math::vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() const -> const std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::positions() -> std::vector<Utils::Math::vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() -> std::vector<Utils::Math::vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() -> std::vector<Utils::Math::vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() -> std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::gatherPositions() const -> std::vector<Utils::Math::vec2d>
{
  std::vector<Utils::Math::vec2d> positions;
  positions.resize(size());
  for (auto [i, pos] : ranges::views::zip(m_globalIndices, m_positions)) {
    positions[i] = pos;
  }
  return positions;
}

auto networkV4::nodes::gatherVelocities() const -> std::vector<Utils::Math::vec2d>
{
  std::vector<Utils::Math::vec2d> velocities;
  velocities.resize(size());
  for (auto [i, vel] : ranges::views::zip(m_globalIndices, m_velocities)) {
    velocities[i] = vel;
  }
  return velocities;
}

auto networkV4::nodes::gatherForces() const -> std::vector<Utils::Math::vec2d>
{
  std::vector<Utils::Math::vec2d> forces;
  forces.resize(size());
  for (auto [i, force] : ranges::views::zip(m_globalIndices, m_forces)) {
    forces[i] = force;
  }
  return forces;
}

auto networkV4::nodes::getNodeMap() const -> const NodeMap
{
  NodeMap nodeMap;
  for (auto [i, index] : ranges::views::enumerate(m_globalIndices)) {
    nodeMap.insert({index, i});
  }
  return nodeMap;
}

// auto networkV4::nodes::nextIndex() -> size_t
//{
//   while (hasIndex(m_nextIndex)) {
//     m_nextIndex++;
//   }
//   return m_nextIndex;
// }

auto networkV4::nodes::hasIndex(size_t _index) const -> bool
{
  return std::find(m_globalIndices.begin(), m_globalIndices.end(), _index)
      != m_globalIndices.end();
}

void networkV4::nodes::zeroVelocity()
{
  std::fill(m_velocities.begin(), m_velocities.end(), Utils::Math::vec2d({0.0, 0.0}));
}

void networkV4::nodes::zeroForce()
{
  // std::fill(m_forces.begin(), m_forces.end(), Utils::Math::vec2d(0.0, 0.0));
  auto view = Utils::spanView(m_forces);
  std::fill(view.begin(), view.end(), 0.0);
}

void networkV4::nodes::pushNode(size_t _globalIndex,
                                const Utils::Math::vec2d& _position,
                                const Utils::Math::vec2d& _velocity,
                                double _mass,
                                const Utils::Math::vec2d& _force)
{
  m_globalIndices.push_back(_globalIndex);
  m_positions.push_back(_position);
  m_velocities.push_back(_velocity);
  m_masses.push_back(_mass);
  m_forces.push_back(_force);
}