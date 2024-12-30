#pragma once

#include <cstdint>
#include <memory>
#include <unordered_map>
#include <variant>
#include <vector>

#include "Core/BreakTypes/None.hpp"
#include "Core/BreakTypes/Strain.hpp"
#include "Core/Forces/Harmonic.hpp"
#include "Core/Forces/Virtual.hpp"
#include "Core/Nodes.hpp"

namespace networkV4
{
namespace bonded
{

struct BondInfo
{
  BondInfo()
      : src {0}
      , dst {0}
  {
  }
  BondInfo(const std::size_t _src, const std::size_t _dst)
      : src {_src}
      , dst {_dst}
  {
  }
  std::size_t src;
  std::size_t dst;
};

struct LocalBondInfo : public BondInfo
{
  LocalBondInfo() = delete;
  LocalBondInfo(const std::size_t _src,
                const std::size_t _dst,
                const std::size_t _index)
      : BondInfo {_src, _dst}
      , index {_index}
  {
  }
  std::size_t index;
};

using bondTypes = std::variant<Forces::VirtualBond, Forces::HarmonicBond>;
using breakTypes = std::variant<BreakTypes::None, BreakTypes::StrainBreak>;

class bondMap
{
public:
  bondMap();
  bondMap(const std::size_t _size);

  void clear();
  void resize(std::size_t _size);
  auto size() const -> const std::size_t;

public:
  void setGlobalInfo(std::size_t _index, size_t _src, size_t _dst);
  void setType(std::size_t _index, const bondTypes& _bond);
  void setBreak(std::size_t _index, const breakTypes& _break);

  auto getGlobalInfo(std::size_t _index) const -> const BondInfo&;
  auto getType(std::size_t _index) const -> const bondTypes&;
  auto getBreak(std::size_t _index) const -> const breakTypes&;

public:
    auto getGlobalInfo() const -> const std::vector<BondInfo>&;
    auto getTypes() const -> const std::vector<bondTypes>&;
    auto getBreaks() const -> const std::vector<breakTypes>&;

public:
  void addTag(std::size_t _index, std::size_t _tag);
  void removeTag(std::size_t _index, std::size_t _tag);
  auto getTagIds(std::size_t _index) const -> const std::vector<std::size_t>&;
  auto hasTag(std::size_t _index, std::size_t _tag) const -> bool const;

private:
  void checkIndex(std::size_t _index) const;
  void checkIndexSet(std::size_t _index) const;

private:
  std::vector<BondInfo> m_globalMap;
  std::vector<bondTypes> m_types;
  std::vector<breakTypes> m_breakTypes;
  std::vector<std::vector<size_t>> m_tags;
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
  auto getLocalInfo(std::size_t _index) const -> const LocalBondInfo&;
  auto getLocal() const -> const std::vector<LocalBondInfo>&;

public:
  auto getMap() const -> const bondMap&;
  auto getMap() -> bondMap&;

public:
  void remap(const NodeMap& _nodeMap);
  void reorder(
      const auto& _order,
      auto fn = [](const auto& _a, const auto& _b) { return _a < _b; });
  void flip();

private:
  void checkIndex(std::size_t _index) const;

private:
  std::vector<LocalBondInfo> m_localMap;

  std::vector<size_t> m_id2Local;

  bondMap m_bonds;
};

}  // namespace bonded
}  // namespace networkV4
