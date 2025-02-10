#pragma once

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Core/BreakTypes/BondedBreak.hpp"
#include "Core/Forces/BondedForces.hpp"
#include "Core/Nodes.hpp"
#include "Misc/Tags/TagStorage.hpp"

namespace networkV4
{
namespace bonded
{

struct BondInfo
{
  BondInfo() = delete;
  BondInfo(const std::size_t _src, const std::size_t _dst)
      : src {_src}
      , dst {_dst}
  {
  }
  std::size_t src;
  std::size_t dst;
};

class bonds
{
public:
  bonds(const size_t _size);

public:
  void clear();
  void reserve(size_t _size);
  size_t size() const;

public:
  auto addBond(size_t _src,
               size_t _dst,
               const bondTypes& _bond = Forces::VirtualBond {},
               const breakTypes& _break = BreakTypes::None {},
               const Utils::Tags::tagFlags& _tags = {}) -> size_t;

public:
  auto getLocalIndex() const -> const std::vector<BondInfo>&;
  auto getGlobalIndex() const -> const std::vector<BondInfo>&;
  auto getTypes() const -> const std::vector<bondTypes>&;
  auto getBreaks() const -> const std::vector<breakTypes>&;
  auto getTags() const -> const Utils::Tags::tagStorage&;

  auto getLocalIndex() -> std::vector<BondInfo>&;
  auto getGlobalIndex() -> std::vector<BondInfo>&;
  auto getTypes() -> std::vector<bondTypes>&;
  auto getBreaks() -> std::vector<breakTypes>&;
  auto getTags() -> Utils::Tags::tagStorage&;

public:
  auto gatherLocalIndex() const -> std::vector<BondInfo> const;
  auto gatherGlobalIndex() const -> std::vector<BondInfo> const;
  auto gatherTypes() const -> std::vector<bondTypes> const;
  auto gatherBreaks() const -> std::vector<breakTypes> const;
  auto gatherTags() const -> Utils::Tags::tagStorage const;

public:
  void remap(const NodeMap& _nodeMap);
  void flipSrcDst();

public:
  template<typename Order>
  void reorder(
      Order _order,
      auto fn = [](const auto& _a, const auto& _b) { return _a < _b; })
  {
    if (_order.size() != size()) {
      throw("bonds::reorder: order size does not match bond size");
    }

    ranges::sort(
        ranges::view::zip(_order, m_localIndex, m_globalIndex, m_indices, m_types, m_breakTypes, m_tags),
        [fn](const auto& _a, const auto& _b)
        { return fn(std::get<0>(_a), std::get<0>(_b)); });
  }

private:
  void boundsCheck(std::size_t _index) const;

private:
  std::vector<BondInfo> m_localIndex;
  std::vector<BondInfo> m_globalIndex;
  std::vector<size_t> m_indices;
  std::vector<bondTypes> m_types;
  std::vector<breakTypes> m_breakTypes;
  Utils::Tags::tagStorage m_tags;
};

}  // namespace bonded
}  // namespace networkV4
