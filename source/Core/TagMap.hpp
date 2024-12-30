#pragma once

#include <string>
#include <unordered_map>

class tagMap
{
public:
  tagMap() {}

public:
  void add(std::size_t _id, const std::string& _name)
  {
    if (m_tagId2Name.find(_id) != m_tagId2Name.end()) {
      throw("tagMap::add: id already exists");
    }

    m_tagId2Name[_id] = _name;
    m_tagName2Id[_name] = _id;
  }

  auto add(const std::string& _name) -> size_t const
  {
    add(updateNextIndex(), _name);
    return m_nextId;
  }

  auto get(const std::string& _name) const -> size_t const
  {
    if (m_tagName2Id.find(_name) == m_tagName2Id.end()) {
      throw("tagMap::get: name does not exist");
    }

    return m_tagName2Id.at(_name);
  }

  auto get(std::size_t _id) const -> const std::string&
  {
    if (m_tagId2Name.find(_id) == m_tagId2Name.end()) {
      throw("tagMap::get: id does not exist");
    }

    return m_tagId2Name.at(_id);
  }

  auto has(const std::string& _name) const -> bool const
  {
    return m_tagName2Id.find(_name) != m_tagName2Id.end();
  }

  auto has(std::size_t _id) const -> bool const
  {
    return m_tagId2Name.find(_id) != m_tagId2Name.end();
  }

private:
  auto updateNextIndex() -> size_t const
  {
    while (m_tagId2Name.find(m_nextId) != m_tagId2Name.end()) {
      m_nextId++;
    }
    return m_nextId;
  }

private:
  std::unordered_map<std::size_t, std::string> m_tagId2Name;
  std::unordered_map<std::string, std::size_t> m_tagName2Id;
  size_t m_nextId;
};