#pragma once

#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>

#include "Misc/Config.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{

using NodeMap = std::unordered_map<size_t, size_t>;

class nodes
{
public:
  nodes(const size_t _size);

public:
  void clear();
  void reserve(size_t _size);
  auto size() const -> size_t;

  void addNode(size_t _globalIndex,
               const Utils::vec2d& _position,
               const Utils::vec2d& _velocity = Utils::vec2d(0.0, 0.0),
               double _mass = 1.0,
               const Utils::vec2d& _force = Utils::vec2d(0.0, 0.0));
  auto addNode(const Utils::vec2d& _position,
               const Utils::vec2d& _velocity = Utils::vec2d(0.0, 0.0),
               double _mass = 1.0) -> size_t;

public:
  auto indices() const -> const std::vector<size_t>&;
  auto positions() const -> const std::vector<Utils::vec2d>&;
  auto velocities() const -> const std::vector<Utils::vec2d>&;
  auto forces() const -> const std::vector<Utils::vec2d>&;
  auto masses() const -> const std::vector<double>&;

  auto positions() -> std::vector<Utils::vec2d>&;
  auto velocities() -> std::vector<Utils::vec2d>&;
  auto forces() -> std::vector<Utils::vec2d>&;
  auto masses() -> std::vector<double>&;

public:
  auto gatherPositions() const -> std::vector<Utils::vec2d>;
  auto getNodeMap() const -> const NodeMap;

public:
  void reorder(
      const auto& _order,
      auto fn = [](const auto& _a, const auto& _b) { return _a < _b; });

public:
  auto nextIndex() -> size_t;
  auto hasIndex(size_t _index) const -> bool;

public:
  void zeroVelocity();
  void zeroForce();

private:
  void pushNode(size_t _globalIndex,
                const Utils::vec2d& _position,
                const Utils::vec2d& _velocity,
                double _mass,
                const Utils::vec2d& _force);

  void checkIndex(std::size_t _index) const;

private:
  std::vector<size_t> m_globalIndices;

  std::vector<Utils::vec2d> m_positions;
  std::vector<Utils::vec2d> m_velocities;
  std::vector<Utils::vec2d> m_forces;

  std::vector<double> m_masses;

  size_t m_nextIndex = 0;
};
}  // namespace networkV4