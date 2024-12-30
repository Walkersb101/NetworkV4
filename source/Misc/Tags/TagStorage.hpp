#pragma once

#include <algorithm>
#include <string>
#include <vector>

#include "Misc/Tags/TagMap.hpp"

namespace Utils
{
namespace Tags
{

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
  void addTag(std::size_t _index, std::size_t _tag)
  {
    boundsCheck(_index);
    if (!contains(m_tags[_index], _tag)) {
      m_tags[_index].push_back(_tag);
    }
  }

  void removeTag(std::size_t _index, std::size_t _tag)
  {
    boundsCheck(_index);
    if (contains(m_tags[_index], _tag)) {
      m_tags[_index].erase(
          std::remove(m_tags[_index].begin(), m_tags[_index].end(), _tag),
          m_tags[_index].end());
    }
  }

  auto getTagIds(std::size_t _index) const -> const std::vector<std::size_t>&
  {
    boundsCheck(_index);
    return m_tags[_index];
  }

  auto getTagStrings(std::size_t _index,
                    const tagMap& _tagMap) const -> std::vector<std::string>
  {
    boundsCheck(_index);
    std::vector<std::string> tagStrings;
    tagStrings.reserve(m_tags[_index].size());
    for (const auto& tag : m_tags[_index]) {
      tagStrings.push_back(_tagMap.get(tag));
    }
    return tagStrings;
  }

  auto hasTag(std::size_t _index, std::size_t _tag) const -> bool const
  {
    boundsCheck(_index);
    return contains(m_tags[_index], _tag);
  }

private:
  void boundsCheck(std::size_t _index) const
  {
    if (_index >= m_tags.size()) {
      throw("tagStorage::boundsCheck: index out of bounds");
    }
  }

  bool contains(const std::vector<size_t>& _tags, std::size_t _tag) const
  {
    return std::find(_tags.begin(), _tags.end(), _tag) != _tags.end();
  }

private:
  std::vector<std::vector<size_t>> m_tags;
};

}  // namespace Tags
}  // namespace Utils