#pragma once

#include <cstdint>
#include <string>
#include <unordered_map>
#include <vector>

#include <range/v3/view/concat.hpp>
#include <range/v3/view/unique.hpp>

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

  auto getTags() const -> std::vector<size_t>
  {
    std::vector<size_t> tags;
    tags.reserve(m_taggedStresses.size());
    for (const auto& [tag, stress] : m_taggedStresses) {
      tags.push_back(tag);
    }
    return tags;
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

public:
  auto has(size_t _tag) const -> bool const
  {
    return m_taggedStresses.find(_tag) != m_taggedStresses.end();
  }

private:
  std::unordered_map<size_t, Utils::tensor2d> m_taggedStresses;
  Utils::tensor2d m_totalStress;
};

inline void merge(stresses& _s1, const stresses& _s2)
{
  _s1.total() += _s2.total();
  auto s1Tags = _s1.getTags();
  auto s2Tags = _s2.getTags();
  for (auto tag : ranges::views::concat(s1Tags, s2Tags) | ranges::view::unique)
  {
    if (_s1.has(tag) && _s2.has(tag)) {
      _s1.get(tag) += _s2.get(tag);
    } else if (_s2.has(tag)) {
      _s1.init(tag);
      _s1.get(tag) = _s2.get(tag);
    }
  }
}

#pragma omp declare reduction(+ : stresses : merge(omp_out, omp_in)) \
    initializer(omp_priv = omp_orig)

}  // namespace networkV4