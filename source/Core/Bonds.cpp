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
  m_bonds.clear();
  m_types.clear();
  m_breakTypes.clear();
    m_tags.clear();
}

void networkV4::bonded::bonds::reserve(std::size_t _size)
{
  m_bonds.reserve(_size);
  m_types.reserve(_size);
  m_breakTypes.reserve(_size);
    m_tags.reserve(_size);
}

auto networkV4::bonded::bonds::size() const -> std::size_t
{
  if (m_bonds.size() != m_types.size() || m_bonds.size() != m_breakTypes.size() || m_bonds.size() != m_tags.size())
  {
    throw("bonds::size: inconsistent sizes");
  }
  return m_bonds.size();
}

auto networkV4::bonded::bonds::addBond(size_t _src,
                                       size_t _dst,
                                       const bondTypes& _bond,
                                       const breakTypes& _break,
                                       const Utils::Tags::tagFlags& _tags) -> size_t
{
  if (_src == _dst) {
    throw("bonds::addBond: src and dst are the same");
  }
  const std::size_t index = m_bonds.size();

  m_bonds.emplace_back(_src, _dst, index);
  m_types.emplace_back(_bond);
  m_breakTypes.emplace_back(_break);
  m_tags.emplace_back(_tags);

  return index;
}

auto networkV4::bonded::bonds::getBonds() const -> const std::vector<BondInfo>&
{
  return m_bonds;
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

auto networkV4::bonded::bonds::getBonds() -> std::vector<BondInfo>&
{
  return m_bonds;
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

auto networkV4::bonded::bonds::gatherBonds() const
    -> std::vector<BondInfo> const
{
  std::vector<BondInfo> bonds;
  bonds.resize(size(), BondInfo(0, 0, 0));
  for (const auto& bond : m_bonds) {
    bonds[bond.index] = bond;
  }
  return bonds;
}

auto networkV4::bonded::bonds::gatherTypes() const
    -> std::vector<bondTypes> const
{
  std::vector<bondTypes> types;
  types.resize(size());
  for (const auto& [type, bond] : ranges::views::zip(m_types, m_bonds)) {
    types[bond.index] = type;
  }
  return types;
}

auto networkV4::bonded::bonds::gatherBreaks() const
    -> std::vector<breakTypes> const
{
  std::vector<breakTypes> breaks;
  breaks.resize(size());
  for (const auto& [brk, bond] : ranges::views::zip(m_breakTypes, m_bonds)) {
    breaks[bond.index] = brk;
  }
  return breaks;
}

auto networkV4::bonded::bonds::gatherTags() const
    -> Utils::Tags::tagStorage const
{
  Utils::Tags::tagStorage tags;
  tags.resize(size());
  for (const auto& [tag, bond] : ranges::views::zip(m_tags, m_bonds)) {
    tags[bond.index] = tag;
  }
  return tags;
}

void networkV4::bonded::bonds::remap(const NodeMap& _nodeMap)
{
  for (auto [i, bond] : ranges::views::enumerate(m_bonds)) {
    m_bonds[i] = BondInfo(_nodeMap.at(bond.src), _nodeMap.at(bond.dst), i);
  }
}

void networkV4::bonded::bonds::flipSrcDst()
{
  for (auto& bond : m_bonds) {
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