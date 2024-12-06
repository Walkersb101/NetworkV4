#pragma once

#include <vector>
#include <cstdint>

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"
#include "Core/box.hpp"
#include <libmorton/morton.h>

namespace networkV4
{
namespace partition
{

class Partition
{
public:
  Partition()
  {
#if defined(_OPENMP)
    const int maxThreads = omp_get_max_threads();
    if (maxThreads > 1) {
      m_partitionsCount = 2 * maxThreads;
    }
#endif
  }

  void assign(const std::vector<Utils::vec2d>& _positions,
              const box& _domain)
  {
    m_partition.clear();
    m_partition.reserve(_positions.size());
    for (const auto& pos : _positions) {
        const auto frac = _domain.toFractional(pos);
        const auto p = static_cast<size_t>(frac.x * m_partitionsCount);
        m_partition.push_back(p);

        const auto x = static_cast<size_t>(frac.x * mortonRes);
        const auto y = static_cast<size_t>(frac.y * mortonRes);
        m_mortonHash.push_back(morton2D(x, y));
    }
  }

private:
  std::size_t m_partitionsCount = 1;

  std::vector<size_t> m_partition;
  std::vector<std::uint_fast64_t> m_mortonHash;
  const size_t mortonRes = 1024;
};

}  // namespace partition
}  // namespace networkV4