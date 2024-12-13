#pragma once

#include <cstdint>
#include <vector>

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <ranges>

#include <libmorton/morton.h>

#include "Core/box.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace partition {

class Partition
{
public:
    Partition(size_t _index, size_t _nodesStart, size_t _nodesEnd, size_t _bondsStart, size_t _bondsEnd)
        : m_index(_index)
        , m_nodesStart(_nodesStart)
        , m_nodesEnd(_nodesEnd)
        , m_bondsStart(_bondsStart)
        , m_bondsEnd(_bondsEnd)
    {
    }
public:
    inline auto index() const -> size_t { return m_index; }

    inline auto nodeStart() const -> size_t { return m_nodesStart; }
    inline auto nodeEnd() const -> size_t { return m_nodesEnd; }
    inline auto nodeCount() const -> size_t { return m_nodesEnd - m_nodesStart; }

    inline auto bondStart() const -> size_t { return m_bondsStart; }
    inline auto bondEnd() const -> size_t { return m_bondsEnd; }
    inline auto bondCount() const -> size_t { return m_bondsEnd - m_bondsStart; }

private:
    size_t m_index;
    size_t m_nodesStart;
    size_t m_nodesEnd;

    size_t m_bondsStart;
    size_t m_bondsEnd;
};

class PartitionGenerator
{
public:
  PartitionGenerator()
  {
#if defined(_OPENMP)
    const int maxThreads = omp_get_max_threads();
    if (maxThreads > 1) {
      m_partitionsCount = 2 * maxThreads;
    }
#endif
  }

  void assign(const std::vector<Utils::vec2d>& _positions, const box& _domain)
  {
    m_partition.clear();
    m_partition.reserve(_positions.size());
    m_mortonHash.clear();
    m_mortonHash.reserve(_positions.size());
    for (const auto& pos : _positions) {
      const auto frac = _domain.toFractional(pos);
      const auto p = static_cast<size_t>(frac.x * m_partitionsCount);
      m_partition.push_back(p);

      const auto x = static_cast<uint_fast32_t>(frac.x * mortonRes);
      const auto y = static_cast<uint_fast32_t>(frac.y * mortonRes);
      const auto hash = libmorton::morton2D_64_encode(x, y);
      m_mortonHash.push_back(hash);
    }
  }

  void sortGlobal();

private:
  std::size_t m_partitionsCount = 1;

  std::vector<size_t> m_partition;
  std::vector<std::uint_fast64_t> m_mortonHash;
  const size_t mortonRes = 1024;
};

}  // namespace partition
}  // namespace networkV4