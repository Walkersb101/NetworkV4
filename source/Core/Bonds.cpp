#include <algorithm>
#include <cmath>
#include <cstddef>
#include <memory>
#include <vector>

#include "Bonds.hpp"

networkV4::bonded::bondMap::bondMap() {}

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

void networkV4::bonded::bondMap::addTag(std::size_t _index, std::size_t _tag)
{
  checkIndex(_index);
  if (!m_tagMap.has(_tag)) {
    throw("bondMap::hasTag: tag does not exist");
  }
  if (!hasTag(_index, _tag)) {
    m_tags[_index].push_back(_tag);
  }
}

void networkV4::bonded::bondMap::addTag(std::size_t _index,
                                        const std::string& _tag)
{
  checkIndex(_index);
  if (!m_tagMap.has(_tag)) {
    throw("bondMap::hasTag: tag does not exist");
  }
  addTag(_index, m_tagMap.get(_tag));
}

void networkV4::bonded::bondMap::removeTag(std::size_t _index, std::size_t _tag)
{
  checkIndex(_index);
  auto& tags = m_tags[_index];
  tags.erase(std::remove(tags.begin(), tags.end(), _tag), tags.end());
}

void networkV4::bonded::bondMap::removeTag(std::size_t _index,
                                           const std::string& _tag)
{
  checkIndex(_index);
  removeTag(_index, m_tagMap.get(_tag));
}

auto networkV4::bonded::bondMap::hasTag(std::size_t _index,
                                        std::size_t _tag) const -> bool const
{
  checkIndex(_index);
  if (!m_tagMap.has(_tag)) {
    throw("bondMap::hasTag: tag does not exist");
  }
  const auto& tags = m_tags[_index];
  return std::find(tags.begin(), tags.end(), _tag) != tags.end();
}

auto networkV4::bonded::bondMap::hasTag(
    std::size_t _index, const std::string& _tag) const -> bool const
{
  checkIndex(_index);
  if (!m_tagMap.has(_tag)) {
    throw("bondMap::hasTag: tag does not exist");
  }
  return hasTag(_index, m_tagMap.get(_tag));
}

void networkV4::bonded::bondMap::checkIndex(std::size_t _index) const
{
  if (_index >= m_types.size()) {
    throw("bondMap::checkIndex: index out of bounds");
  }
}