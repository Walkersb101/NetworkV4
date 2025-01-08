#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "Misc/Tags/TagMap.hpp"

namespace Utils
{
namespace Tags
{

using tagFlags = std::bitset<NUM_TAGS>;
using tagStorage = std::vector<tagFlags>;

inline auto hasTag(const tagFlags& _tags, const tagFlags& _tag) -> bool
{
  return (_tags & _tag) == _tag;
}

/*
class tagStorage
{
public:
  tagStorage() {}
  tagStorage(const std::size_t _size)
      : m_tags(_size)
  {
  }

public:
  void clear() { m_tags.clear(); }

  void resize(std::size_t _size)
  {
    m_tags.resize(_size);
    // TODO: add log if resizing to smaller size
  }

  auto size() const -> const std::size_t { return m_tags.size(); }

public:
  void addTags(std::size_t _index, const std::bitset<NUM_TAGS>& _tags)
  {
    boundsCheck(_index);
    m_tags.at(_index) |= _tags;
  }

  void removeTag(std::size_t _index, const std::bitset<NUM_TAGS>& _tags)
  {
    boundsCheck(_index);
    m_tags.at(_index) &= ~_tags;
  }

  auto getTagIds(std::size_t _index) const -> const std::bitset<NUM_TAGS>&
  {
    boundsCheck(_index);
    return m_tags.at(_index);
  }

  auto hasTag(std::size_t _index, std::size_t _tag) const -> bool
  {
    boundsCheck(_index);
    tagCheck(_tag);
    return m_tags.at(_index).test(_tag);
  }

  auto data() -> std::vector<std::bitset<NUM_TAGS>>& { return m_tags; }
    auto data() const -> const std::vector<std::bitset<NUM_TAGS>>& { return
m_tags; }

private:
  void boundsCheck(std::size_t _index) const
  {
    if (_index >= size()) {
      throw("bonds::checkIndex: index out of bounds");  // TODO: should this
                                                        // throw?
    }
  }

  void tagCheck(std::size_t _tag) const
  {
    if (_tag >= NUM_TAGS) {
      throw("bonds::checkTag: tag out of bounds");  // TODO: should this
                                                    // throw?
    }
  }

private:
  std::vector<std::bitset<NUM_TAGS>> m_tags;
};
*/

}  // namespace Tags
}  // namespace Utils