#pragma once

#include <cstdint>
#include <vector>

#include "Misc/Math/Tensor2.hpp"

namespace networkV4
{

class stresses
{
public:
  stresses() {}

public:
  void clear()
  {
    m_stored.reset();
    zero();
  }

  void init(const std::bitset<NUM_TAGS>& _tags,
            const Utils::tensor2d& _value = Utils::tensor2d())
  {
    std::bitset<NUM_TAGS> diff = _tags & ~m_stored;
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      // TODO : ADD LOG if already set
      if (diff.test(i)) {
        init(i, _value);
      }
    }
  }

  auto getFirst(const std::bitset<NUM_TAGS>& _tags) -> Utils::tensor2d&
  {
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      if (_tags.test(i)) {
        return m_values.at(i);
      }
    }
    throw std::invalid_argument("Tag not found");
  }

  auto getFirst(const std::bitset<NUM_TAGS>& _tags) const
      -> const Utils::tensor2d&
  {
    return const_cast<stresses*>(this)->getFirst(_tags);
  }

  auto getall(const std::bitset<NUM_TAGS>& _tags)
      -> std::vector<Utils::tensor2d>
  {
    std::vector<Utils::tensor2d> values;
    values.reserve(_tags.count());
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      if (_tags.test(i)) {
        values.emplace_back(m_values.at(i));
      }
    }
    return values;
  }

  auto getInitilised() -> const std::bitset<NUM_TAGS>& { return m_stored; }

  auto total() -> Utils::tensor2d& { return m_totalStress; }
  auto total() const -> const Utils::tensor2d& { return m_totalStress; }

  void zero()
  {
    std::fill(m_values.begin(), m_values.end(), Utils::tensor2d());
    m_totalStress = Utils::tensor2d();
  }

  void distribute(const Utils::tensor2d& _stress,
                  const std::bitset<NUM_TAGS>& _tags)
  {
    std::bitset<NUM_TAGS> overlap = m_stored & _tags;
    m_totalStress += _stress;
    for (size_t i = 0; i < NUM_TAGS; ++i) {
      if (overlap.test(i)) {
        m_values.at(i) += _stress;
      }
    }
  }

private:
  void init(size_t _index, const Utils::tensor2d& _value)
  {
    m_stored.set(_index);
    m_values.at(_index) = _value;
  }

private:
  std::bitset<NUM_TAGS> m_stored;
  std::array<Utils::tensor2d, NUM_TAGS> m_values;
  Utils::tensor2d m_totalStress;

private:
  friend void merge(stresses& _s1, const stresses& _s2);
};

inline void merge(stresses& _s1, const stresses& _s2)
{
  _s1.m_totalStress += _s2.m_totalStress;

  auto s1init = _s1.m_stored;
  auto s2init = _s2.m_stored;

  for (size_t i = 0; i < NUM_TAGS; ++i) {
    if (s2init.test(i)) {
      if (s1init.test(i)) {
        _s1.m_values[i] += _s2.m_values[i];
      } else {
        _s1.m_stored.set(i);
        _s1.m_values[i] = _s2.m_values[i];
      }
    }
  }
}

#pragma omp declare reduction(+ : stresses : merge(omp_out, omp_in)) \
    initializer(omp_priv = omp_orig)

}  // namespace networkV4