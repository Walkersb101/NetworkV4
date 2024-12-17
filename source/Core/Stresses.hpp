#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include "Misc/Tensor2.hpp"

namespace networkV4
{

class stresses
{
public:
  stresses() {}

public:
  void clear() { m_taggedStresses.clear(); }

  void init(size_t _tag)
  {
    if (!has(_tag)) {
      m_taggedStresses[_tag] = Utils::tensor2d(0.0, 0.0, 0.0, 0.0);
    }  // TODO: else throw? Or Log?
  }
  auto get(size_t _tag) -> Utils::tensor2d&
  {
    return m_taggedStresses.at(_tag);
  }
  auto get(size_t _tag) const -> const Utils::tensor2d&
  {
    return m_taggedStresses.at(_tag);
  }

  auto total() -> Utils::tensor2d& { return m_totalStress; }
  auto total() const -> const Utils::tensor2d& { return m_totalStress; }

  void zero()
  {
    for (auto& [tag, stress] : m_taggedStresses) {
      stress.set(0.0, 0.0, 0.0, 0.0);
    }
    m_totalStress.set(0.0, 0.0, 0.0, 0.0);
  }

  void distribute(const Utils::tensor2d& _stress,
                  const std::vector<std::size_t>& _tags)
  {
    m_totalStress += _stress;
    for (const auto& tag : _tags) {
      if (has(tag)) {
        m_taggedStresses[tag] += _stress;
      }
    }
  }

private:
  auto has(size_t _tag) const -> bool const
  {
    return m_taggedStresses.find(_tag) != m_taggedStresses.end();
  }

private:
  std::unordered_map<size_t, Utils::tensor2d> m_taggedStresses;
  Utils::tensor2d m_totalStress;
};
}  // namespace networkV4