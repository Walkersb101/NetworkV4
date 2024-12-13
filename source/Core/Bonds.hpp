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
#include "TagMap.hpp"

namespace networkV4
{
namespace bonded
{

struct BondInfo
{
  const std::size_t src;
  const std::size_t dst;
  const std::size_t id;
};

class bonds
{

public:


private:
    std::vector<BondInfo> m_globalMap;
    std::vector<BondInfo> m_localMap;

    std::vector<std::size_t> m_local2Global;

    bondMap m_bonds;
};

using bondTypes = std::variant<Forces::VirtualBond, Forces::HarmonicBond>;
using breakTypes = std::variant<BreakTypes::None, BreakTypes::StrainBreak>;

class bondMap
{
public:
  bondMap();

  void clear();
  void resize(std::size_t _size);

  void set(std::size_t _index, const bondTypes& _bond);
  void set(std::size_t _index, const breakTypes& _break);
  void remove(std::size_t _index);

  auto size() const -> std::size_t const;

  void addTag(std::size_t _index, std::size_t _tag);
  void addTag(std::size_t _index, const std::string& _tag);
  void removeTag(std::size_t _index, std::size_t _tag);
  void removeTag(std::size_t _index, const std::string& _tag);

  auto hasTag(std::size_t _index, std::size_t _tag) const -> bool const;
  auto hasTag(std::size_t _index, const std::string& _tag) const -> bool const;

private:
  void checkIndex(std::size_t _index) const;

private:
  std::vector<bondTypes> m_types;
  std::vector<breakTypes> m_breakTypes;
  std::vector<std::vector<size_t>> m_tags;

  tagMap m_tagMap;
};

}  // namespace bonded
}  // namespace networkV4
