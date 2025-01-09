#pragma once

#include <string>
#include <vector>

namespace Utils
{
namespace Tags
{

#define BROKEN_TAG_INDEX 0

class tagMap
{
public:
  tagMap() 
  {
  }

public:
  auto size() const -> size_t const { return NUM_TAGS; }

  auto add(const std::string& _name) -> std::bitset<NUM_TAGS>
  {
    size_t index = firstEmpty();
    m_tagNames[index] = _name;
    m_setTags.set(index);
    return index2Bitset(index);
  }

  auto get(const std::string& _name) const -> std::bitset<NUM_TAGS>
  {
    auto it = std::find(m_tagNames.begin(), m_tagNames.end(), _name);
    if (it == m_tagNames.end()) {
      throw "tagMap::get: name does not exist";
    }
    auto index = std::distance(m_tagNames.begin(), it);
    return std::bitset<NUM_TAGS>(1 << index);
  }

  auto get(std::size_t _id) const -> const std::string&
  {
    checkTag(_id);
    return m_tagNames[_id];
  }

  auto get(const std::bitset<NUM_TAGS>& _tags) const -> std::vector<std::string>
  {
    std::vector<std::string> names;
    names.reserve(_tags.count());
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      if (_tags.test(i)) {
        names.emplace_back(m_tagNames[i]);
      }
    }
    return names;
  }

  auto has(const std::string& _name) const -> bool const
  {
    return std::find(m_tagNames.begin(), m_tagNames.end(), _name)
        != m_tagNames.end();
  }

  auto has(std::size_t _id) const -> bool const
  {
    checkTag(_id);
    return m_setTags.test(_id);
  }

  auto getVecs() const -> const std::array<std::string, NUM_TAGS>&
  {
    return m_tagNames;
  }

private:
  void checkTag(std::size_t _tag) const
  {
    if (_tag >= NUM_TAGS) {
      throw "tagMap::checkTag: tag out of bounds";
    }
  }

  auto firstEmpty() const -> size_t const
  {
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      if (!m_setTags.test(i)) {
        return i;
      }
    }

    throw "tagMap::firstEmpty: TagMap is full";
  }

  auto index2Bitset(std::size_t _index) const -> std::bitset<NUM_TAGS>
  {
    checkTag(_index);
    return std::bitset<NUM_TAGS>(1 << _index);
  }

private:
  std::array<std::string, NUM_TAGS> m_tagNames;
  std::bitset<NUM_TAGS> m_setTags;
};

}  // namespace Tags
}  // namespace Utils