#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "Bonds.hpp"

#include <range/v3/view/enumerate.hpp>
#include <range/v3/algorithm.hpp>


// -------------------------------------------------------------------------- //

networkV4::bonded::bondMap::bondMap() {}

networkV4::bonded::bondMap::bondMap(std::size_t _size)
{
  resize(_size);
}

void networkV4::bonded::bondMap::clear()
{
  m_types.clear();
  m_breakTypes.clear();
  m_tags.clear();
}

void networkV4::bonded::bondMap::resize(std::size_t _size)
{
  m_types.resize(_size);
  m_breakTypes.resize(_size);
  m_tags.resize(_size);
}

const std::size_t networkV4::bonded::bondMap::size() const
{
  if (m_types.size() != m_breakTypes.size() || m_types.size() != m_tags.size())
  {
    throw("bondMap::size: inconsistent sizes");
  }
  return m_types.size();
}

void networkV4::bonded::bondMap::setGlobalInfo(std::size_t _index,
                                               size_t _src,
                                               size_t _dst)
{
  checkIndex(_index);
  m_globalMap[_index] = BondInfo(_src, _dst);
}

void networkV4::bonded::bondMap::setType(std::size_t _index,
                                         const bondTypes& _bond)
{
  checkIndexSet(_index);
  m_types[_index] = _bond;
}

void networkV4::bonded::bondMap::setBreak(std::size_t _index,
                                          const breakTypes& _break)
{
  checkIndexSet(_index);
  m_breakTypes[_index] = _break;
}

auto networkV4::bonded::bondMap::getGlobalInfo(std::size_t _index) const
    -> const BondInfo&
{
  checkIndex(_index);
  return m_globalMap[_index];
}

auto networkV4::bonded::bondMap::getType(std::size_t _index) const
    -> const bondTypes&
{
  checkIndex(_index);
  return m_types[_index];
}

auto networkV4::bonded::bondMap::getBreak(std::size_t _index) const
    -> const breakTypes&
{
  checkIndex(_index);
  return m_breakTypes[_index];
}

void networkV4::bonded::bondMap::addTag(std::size_t _index, std::size_t _tag)
{
  checkIndex(_index);
  if (!hasTag(_index, _tag)) {
    m_tags[_index].push_back(_tag);
  }
}

void networkV4::bonded::bondMap::removeTag(std::size_t _index, std::size_t _tag)
{
  checkIndex(_index);
  if (!hasTag(_index, _tag)) {
    return;  // TODO: add log if trying to remove tag that does not exist
  }
  auto& tags = m_tags[_index];
  tags.erase(std::remove(tags.begin(), tags.end(), _tag), tags.end());
}

auto networkV4::bonded::bondMap::getTagIds(std::size_t _index) const
    -> const std::vector<std::size_t>&
{
  checkIndex(_index);
  return m_tags[_index];
}

auto networkV4::bonded::bondMap::hasTag(std::size_t _index,
                                        std::size_t _tag) const -> bool const
{
  checkIndex(_index);
  const auto& tags = m_tags[_index];
  return std::find(tags.begin(), tags.end(), _tag) != tags.end();
}

void networkV4::bonded::bondMap::checkIndex(std::size_t _index) const
{
  if (_index >= m_types.size()) {
    throw("bondMap::checkIndex: index out of bounds");  // TODO: maybe this
                                                        // shouldn't throw?
  }
}

void networkV4::bonded::bondMap::checkIndexSet(std::size_t _index) const
{
  checkIndex(_index);
  const auto& binfo = m_globalMap[_index];
  if (binfo.src == 0 && binfo.dst == 0) {
    throw("bondMap::checkIndexSet: index not set");
  }
}

// -------------------------------------------------------------------------- //

networkV4::bonded::bonds::bonds(const size_t _size)
{
  reserve(_size);
}

void networkV4::bonded::bonds::clear()
{
  m_localMap.clear();
  m_id2Local.clear();
  m_bonds.clear();
}

void networkV4::bonded::bonds::reserve(std::size_t _size)
{
  m_localMap.reserve(_size);
  m_id2Local.reserve(_size);
  m_bonds.resize(_size);
}

auto networkV4::bonded::bonds::size() const -> std::size_t
{
  if (m_localMap.size() != m_id2Local.size()
      || m_localMap.size() != m_bonds.size())
  {
    throw("bonds::size: inconsistent sizes");
  }
  return m_localMap.size();
}

void networkV4::bonded::bonds::addBond(size_t _src,
                                       size_t _dst,
                                       const bondTypes& _bond,
                                       const breakTypes& _break)
{
  if (_src == _dst) {
    throw("bonds::addBond: src and dst are the same");
  }
  const std::size_t index = m_bonds.size();
  m_localMap.emplace_back(_src, _dst, index);
  m_id2Local.push_back(index);
  if (m_bonds.size() <= index) {
    m_bonds.resize(index + 1);
  }
  m_bonds.setGlobalInfo(index, _src, _dst);
  m_bonds.setType(index, _bond);
  m_bonds.setBreak(index, _break);
}

auto networkV4::bonded::bonds::getLocalInfo(std::size_t _index) const
    -> const LocalBondInfo&
{
  checkIndex(_index);
  return m_localMap[m_id2Local[_index]];
}

auto networkV4::bonded::bonds::getLocal() const
    -> const std::vector<LocalBondInfo>&
{
  return m_localMap;
}

auto networkV4::bonded::bonds::getMap() const
    -> const networkV4::bonded::bondMap&
{
  return m_bonds;
}

auto networkV4::bonded::bonds::getMap() -> networkV4::bonded::bondMap&
{
  return m_bonds;
}

void networkV4::bonded::bonds::remap(const NodeMap& _nodeMap)
{
  for (auto [i, bond] : ranges::views::enumerate(m_localMap)) {
    m_localMap[i] =
        LocalBondInfo(_nodeMap.at(bond.src), _nodeMap.at(bond.dst), i);
  }
}

void networkV4::bonded::bonds::reorder(const auto& _order, auto fn)
{
  if (_order.size() != size()) {
    throw("bonds::reorder: order size does not match bond size");
  }

  ranges::sort(ranges::view::zip(_order, m_localMap),
               [fn](const auto& _a, const auto& _b)
               { return fn(std::get<0>(_a), std::get<0>(_b)); });
}

void networkV4::bonded::bonds::flip()
{
  for (auto& bond : m_localMap) {
    if (bond.src > bond.dst) {
      std::swap(bond.src, bond.dst);
    }
  }
}

void networkV4::bonded::bonds::checkIndex(std::size_t _index) const
{
  if (_index >= m_localMap.size()) {
    throw("bonds::checkIndex: index out of bounds");
  }
}