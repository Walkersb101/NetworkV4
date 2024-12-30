#pragma once

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <vector>

#include "Core/BreakTypes/BondedBreak.hpp"
#include "Core/Forces/BondedForces.hpp"
#include "Core/Nodes.hpp"

namespace networkV4
{
namespace bonded
{

struct BondInfo
{
  BondInfo() = delete;
  BondInfo(const std::size_t _src,
           const std::size_t _dst,
           const std::size_t _index)
      : src {_src}
      , dst {_dst}
      , index {_index}
  {
  }
  std::size_t src;
  std::size_t dst;
  std::size_t index;
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
               const breakTypes& _break = BreakTypes::None {}) -> size_t;

public:
  auto getBonds() const -> const std::vector<BondInfo>&;
  auto getTypes() const -> const std::vector<bondTypes>&;
  auto getBreaks() const -> const std::vector<breakTypes>&;

  auto getBonds() -> std::vector<BondInfo>&;
  auto getTypes() -> std::vector<bondTypes>&;
  auto getBreaks() -> std::vector<breakTypes>&;

public:
  auto gatherBonds() const -> std::vector<BondInfo>;
  auto gatherTypes() const -> std::vector<bondTypes>;
  auto gatherBreaks() const -> std::vector<breakTypes>;

public:
  void remap(const NodeMap& _nodeMap);
  void reorder(
      const auto& _order,
      auto fn = [](const auto& _a, const auto& _b) { return _a < _b; });
  void flipSrcDst();

private:
  void boundsCheck(std::size_t _index) const;

private:
  std::vector<BondInfo> m_bonds;
  std::vector<bondTypes> m_types;
  std::vector<breakTypes> m_breakTypes;
};

}  // namespace bonded
}  // namespace networkV4
