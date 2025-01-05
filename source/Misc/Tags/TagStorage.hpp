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
    m_tags.reserve(_size);
    // TODO: add log if resizing to smaller size
  }

  auto size() const -> const std::size_t { return m_tags.size(); }

public:
  void initIndex(std::size_t _index) { m_tags.insert({_index, {}}); }

  void addTag(std::size_t _index, std::size_t _tag)
  {
    if (keyExists(_index) && !contains(m_tags.at(_index), _tag)) {
      m_tags.at(_index).push_back(_tag);
    } else {
      m_tags.insert({_index, {_tag}});
    }
  }

  void removeTag(std::size_t _index, std::size_t _tag)
  {
    if (keyExists(_index) && contains(m_tags[_index], _tag)) {
      m_tags[_index].erase(
          std::remove(m_tags[_index].begin(), m_tags[_index].end(), _tag),
          m_tags[_index].end());
    }
  }

  auto getTagIds(std::size_t _index) const -> std::vector<std::size_t>
  {
    if (keyExists(_index)) {
      return m_tags.at(_index);
    } else {
        return {};
    }
  }

  auto getTagStrings(std::size_t _index,
                     const tagMap& _tagMap) const -> std::vector<std::string>
  {
    if (!keyExists(_index)) {
      return {};
    }
    std::vector<std::string> tagStrings;
    tagStrings.reserve(m_tags.at(_index).size());
    for (const auto& tag : m_tags.at(_index)) {
      tagStrings.push_back(_tagMap.get(tag));
    }
    return tagStrings;
  }

  auto hasTag(std::size_t _index, std::size_t _tag) const -> bool const
  {
    return keyExists(_index) && contains(m_tags.at(_index), _tag);
  }

private:
  bool keyExists(std::size_t _index) const { return m_tags.contains(_index); }

  bool contains(const std::vector<size_t>& _tags, std::size_t _tag) const
  {
    return std::find(_tags.begin(), _tags.end(), _tag) != _tags.end();
  }

private:
  std::unordered_map<size_t, std::vector<size_t>> m_tags;
};

}  // namespace Tags
}  // namespace Utils