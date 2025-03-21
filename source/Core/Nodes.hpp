#pragma once

#include <cstdint>
#include <set>
#include <unordered_map>
#include <vector>
#include <range/v3/all.hpp>

#include "Misc/Config.hpp"
#include "Misc/Math/Vector.hpp"

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

  //void addNode(size_t _globalIndex,
  //             const Utils::Math::vec2d& _position,
  //             const Utils::Math::vec2d& _velocity = Utils::Math::vec2d(0.0, 0.0),
  //             double _mass = 1.0,
  //             const Utils::Math::vec2d& _force = Utils::Math::vec2d(0.0, 0.0));
  auto addNode(const Utils::Math::vec2d& _position,
               const Utils::Math::vec2d& _velocity = Utils::Math::vec2d({0.0, 0.0}),
               double _mass = 1.0) -> size_t;

public:
  auto indices() const -> const std::vector<size_t>&;
  auto positions() const -> const std::vector<Utils::Math::vec2d>&;
  auto velocities() const -> const std::vector<Utils::Math::vec2d>&;
  auto forces() const -> const std::vector<Utils::Math::vec2d>&;
  auto masses() const -> const std::vector<double>&;

  auto positions() -> std::vector<Utils::Math::vec2d>&;
  auto velocities() -> std::vector<Utils::Math::vec2d>&;
  auto forces() -> std::vector<Utils::Math::vec2d>&;
  auto masses() -> std::vector<double>&;

public:
  auto gatherPositions() const -> std::vector<Utils::Math::vec2d>;
  auto gatherVelocities() const -> std::vector<Utils::Math::vec2d>;
  auto gatherForces() const -> std::vector<Utils::Math::vec2d>;
  auto getNodeMap() const -> const NodeMap;

public:
  template<typename Order>
  void reorder(
      Order _order,
      auto fn = [](const auto& _a, const auto& _b) { return _a < _b; })
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

public:
  //auto nextIndex() -> size_t;
  auto hasIndex(size_t _index) const -> bool;

public:
  void zeroVelocity();
  void zeroForce();

private:
  void pushNode(size_t _globalIndex,
                const Utils::Math::vec2d& _position,
                const Utils::Math::vec2d& _velocity,
                double _mass,
                const Utils::Math::vec2d& _force);

  void checkIndex(std::size_t _index) const;

private:
  std::vector<size_t> m_globalIndices;

  std::vector<Utils::Math::vec2d> m_positions;
  std::vector<Utils::Math::vec2d> m_velocities;
  std::vector<Utils::Math::vec2d> m_forces;

  std::vector<double> m_masses;

  //size_t m_nextIndex = 0;
};
}  // namespace networkV4