#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "Bonds.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view/enumerate.hpp>

// -------------------------------------------------------------------------- //

networkV4::bonded::bonds::bonds(const size_t _size)
{
  reserve(_size);
}

void networkV4::bonded::bonds::clear()
{
  m_localIndex.clear();
  m_globalIndex.clear();
  m_indices.clear();
  m_types.clear();
  m_breakTypes.clear();
  m_tags.clear();
}

void networkV4::bonded::bonds::reserve(std::size_t _size)
{
  m_localIndex.reserve(_size);
  m_globalIndex.reserve(_size);
  m_indices.reserve(_size);
  m_types.reserve(_size);
  m_breakTypes.reserve(_size);
  m_tags.reserve(_size);
}

auto networkV4::bonded::bonds::size() const -> std::size_t
{
  if (m_localIndex.size() != m_types.size()
      || m_localIndex.size() != m_breakTypes.size()
      || m_localIndex.size() != m_tags.size())
  {
    throw("bonds::size: inconsistent sizes");
  }
  return m_indices.size();
}

auto networkV4::bonded::bonds::addBond(size_t _src,
                                       size_t _dst,
                                       const bondTypes& _bond,
                                       const breakTypes& _break,
                                       const Utils::Tags::tagFlags& _tags)
    -> size_t
{
  if (_src == _dst) {
    throw("bonds::addBond: src and dst are the same");
  }
  const std::size_t index = m_indices.size();

  m_localIndex.emplace_back(_src, _dst);
  m_globalIndex.emplace_back(_src, _dst);
  m_indices.emplace_back(index);
  m_types.emplace_back(_bond);
  m_breakTypes.emplace_back(_break);
  m_tags.emplace_back(_tags);

  return index;
}

auto networkV4::bonded::bonds::getLocalIndex() const
    -> const std::vector<BondInfo>&
{
  return m_localIndex;
}

auto networkV4::bonded::bonds::getGlobalIndex() const
    -> const std::vector<BondInfo>&
{
  return m_globalIndex;
}

auto networkV4::bonded::bonds::getTypes() const -> const std::vector<bondTypes>&
{
  return m_types;
}

auto networkV4::bonded::bonds::getBreaks() const
    -> const std::vector<breakTypes>&
{
  return m_breakTypes;
}

auto networkV4::bonded::bonds::getTags() const -> const Utils::Tags::tagStorage&
{
  return m_tags;
}

auto networkV4::bonded::bonds::getLocalIndex() -> std::vector<BondInfo>&
{
  return m_localIndex;
}

auto networkV4::bonded::bonds::getGlobalIndex() -> std::vector<BondInfo>&
{
  return m_globalIndex;
}

auto networkV4::bonded::bonds::getTypes() -> std::vector<bondTypes>&
{
  return m_types;
}

auto networkV4::bonded::bonds::getBreaks() -> std::vector<breakTypes>&
{
  return m_breakTypes;
}

auto networkV4::bonded::bonds::getTags() -> Utils::Tags::tagStorage&
{
  return m_tags;
}

auto networkV4::bonded::bonds::gatherLocalIndex() const
    -> std::vector<BondInfo> const
{
  std::vector<BondInfo> bonds;
  bonds.resize(size(), BondInfo(0, 0));
  for (const auto& [bond, index] : ranges::views::zip(m_localIndex, m_indices))
  {
    bonds[index] = bond;
  }
  return bonds;
}

auto networkV4::bonded::bonds::gatherGlobalIndex() const
    -> std::vector<BondInfo> const
{
  std::vector<BondInfo> bonds;
  bonds.resize(size(), BondInfo(0, 0));
  for (const auto& [bond, index] : ranges::views::zip(m_globalIndex, m_indices))
  {
    bonds[index] = bond;
  }
  return bonds;
}

auto networkV4::bonded::bonds::gatherTypes() const
    -> std::vector<bondTypes> const
{
  std::vector<bondTypes> types;
  types.resize(size());
  for (const auto& [type, index] : ranges::views::zip(m_types, m_indices)) {
    types[index] = type;
  }
  return types;
}

auto networkV4::bonded::bonds::gatherBreaks() const
    -> std::vector<breakTypes> const
{
  std::vector<breakTypes> breaks;
  breaks.resize(size());
  for (const auto& [brk, index] : ranges::views::zip(m_breakTypes, m_indices)) {
    breaks[index] = brk;
  }
  return breaks;
}

auto networkV4::bonded::bonds::gatherTags() const
    -> Utils::Tags::tagStorage const
{
  Utils::Tags::tagStorage tags;
  tags.resize(size());
  for (const auto& [tag, index] : ranges::views::zip(m_tags, m_indices)) {
    tags[index] = tag;
  }
  return tags;
}

void networkV4::bonded::bonds::remap(const NodeMap& _nodeMap)
{
  for (auto [lI, gI] : ranges::views::zip(m_localIndex, m_globalIndex)) {
    lI = BondInfo(_nodeMap.at(gI.src), _nodeMap.at(gI.dst));
  }
}

void networkV4::bonded::bonds::flipSrcDst()
{
  for (auto& bond : m_localIndex) {
    if (bond.src > bond.dst) {
      std::swap(bond.src, bond.dst);
    }
  }
}

void networkV4::bonded::bonds::boundsCheck(std::size_t _index) const
{
  if (_index >= size()) {
    throw("bonds::checkIndex: index out of bounds");
  }
}