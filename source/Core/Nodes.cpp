#include <algorithm>
#include <stdexcept>

#include "Core/Nodes.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view/enumerate.hpp>
#include <range/v3/view/zip.hpp>

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

void networkV4::nodes::addNode(size_t _globalIndex,
                               const Utils::vec2d& _position,
                               const Utils::vec2d& _velocity,
                               double _mass,
                               const Utils::vec2d& _force)
{
  if (hasIndex(_globalIndex)) {
    throw std::runtime_error("Node with index " + std::to_string(_globalIndex)
                             + " already exists");
  }
  // TODO: check if tags have been set
  pushNode(_globalIndex, _position, _velocity, _mass, _force);
}

auto networkV4::nodes::addNode(const Utils::vec2d& _position,
                               const Utils::vec2d& _velocity,
                               double _mass) -> size_t
{
  // TODO: check if tags have been set
  size_t index = nextIndex();
  pushNode(index, _position, _velocity, _mass, Utils::vec2d(0.0, 0.0));
  m_nextIndex++;
  return index;
}

auto networkV4::nodes::indices() const -> const std::vector<size_t>&
{
  return m_globalIndices;
}

auto networkV4::nodes::positions() const -> const std::vector<Utils::vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() const -> const std::vector<Utils::vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() const -> const std::vector<Utils::vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() const -> const std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::positions() -> std::vector<Utils::vec2d>&
{
  return m_positions;
}

auto networkV4::nodes::velocities() -> std::vector<Utils::vec2d>&
{
  return m_velocities;
}

auto networkV4::nodes::forces() -> std::vector<Utils::vec2d>&
{
  return m_forces;
}

auto networkV4::nodes::masses() -> std::vector<double>&
{
  return m_masses;
}

auto networkV4::nodes::gatherPositions() const -> std::vector<Utils::vec2d>
{
  std::vector<Utils::vec2d> positions;
  positions.resize(size());
  for (auto [i, pos] : ranges::views::zip(m_globalIndices, m_positions)) {
    positions[i] = pos;
  }
  return positions;
}

auto networkV4::nodes::getNodeMap() const -> const NodeMap
{
  NodeMap nodeMap;
  for (auto [i, index] : ranges::views::enumerate(m_globalIndices)) {
    nodeMap.insert({index, i});
  }
  return nodeMap;
}

void networkV4::nodes::reorder(const auto& _order, auto fn)
{
  if (_order.size() != size()) {
    throw std::runtime_error(
        "nodes::reorder: order size does not match node size");
  }

  ranges::sort(ranges::view::zip(_order,
                                 m_globalIndices,
                                 m_positions,
                                 m_velocities,
                                 m_forces,
                                 m_masses),
               [fn](const auto& _a, const auto& _b)
               { return fn(std::get<0>(_a), std::get<0>(_b)); });
}

void networkV4::nodes::addTag(std::size_t _index, std::size_t _tag)
{
  if (hasTag(_index, _tag)) {
    return;  // TODO: add log if trying to add tag that already exists
  }
  if (m_tags.find(_index) == m_tags.end()) {
    // tagList for this index does not exist
    if (!hasIndex(_index)) {  // check if node with this index exists
      throw std::runtime_error("nodes::addTag: index does not exist");
    }
    m_tags.insert({_index, {_tag}});
  } else {
    m_tags.at(_index).push_back(_tag);
  }
}

void networkV4::nodes::removeTag(std::size_t _index, std::size_t _tag)
{
  if (!hasTag(_index, _tag)) {
    return;  // TODO: add log if trying to remove tag that does not exist
  }
  auto& tags = m_tags.at(_index);
  tags.erase(std::remove(tags.begin(), tags.end(), _tag), tags.end());
}

auto networkV4::nodes::hasTag(std::size_t _index,
                              std::size_t _tag) const -> bool const
{
  if (m_tags.find(_index) == m_tags.end()) {
    return false;
  }
  const auto& tags = m_tags.at(_index);
  return std::find(tags.begin(), tags.end(), _tag) != tags.end();
}

auto networkV4::nodes::nextIndex() -> size_t
{
  while (hasIndex(m_nextIndex)) {
    m_nextIndex++;
  }
  return m_nextIndex;
}

auto networkV4::nodes::hasIndex(size_t _index) const -> bool
{
  return std::find(m_globalIndices.begin(), m_globalIndices.end(), _index)
      != m_globalIndices.end();
}

void networkV4::nodes::zeroVelocity()
{
  std::fill(m_velocities.begin(), m_velocities.end(), Utils::vec2d(0.0, 0.0));
}

void networkV4::nodes::zeroForce()
{
  std::fill(m_forces.begin(), m_forces.end(), Utils::vec2d(0.0, 0.0));
}

void networkV4::nodes::pushNode(size_t _globalIndex,
                                const Utils::vec2d& _position,
                                const Utils::vec2d& _velocity,
                                double _mass,
                                const Utils::vec2d& _force)
{
  m_globalIndices.push_back(_globalIndex);
  m_positions.push_back(_position);
  m_velocities.push_back(_velocity);
  m_masses.push_back(_mass);
  m_forces.push_back(_force);
}