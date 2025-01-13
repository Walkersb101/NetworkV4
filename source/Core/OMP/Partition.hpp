#pragma once

#include <cstdint>
#include <stdexcept>
#include <vector>

#if defined(_OPENMP)
#  include <omp.h>
#endif

#include <libmorton/morton.h>
#include <range/v3/algorithm.hpp>
#include <range/v3/view.hpp>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Core/box.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Math/Vector.hpp"

namespace networkV4
{
namespace partition
{

class Partition
{
public:
  Partition(size_t _index,
            size_t _nodesStart,
            size_t _nodesEnd,
            size_t _bondsStart,
            size_t _bondsEnd)
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

using Partitions = std::vector<Partition>;

class PartitionGenerator
{
public:
  PartitionGenerator()
  {
#if defined(_OPENMP)
    int maxThreads = omp_get_max_threads();
    // if (maxThreads > 1) {
    //maxThreads = 2;
    m_partitionsCount = 2 * maxThreads;
    m_passes = 2;
    //}
#endif
  }

  void assignNodes(const std::vector<Utils::Math::vec2d>& _positions,
                   const box& _domain)
  {
    m_partition.clear();
    m_partition.reserve(_positions.size());
    m_mortonHash.clear();
    m_mortonHash.reserve(_positions.size());

    const auto partitionSize = 1.0 / m_partitionsCount;
    for (const auto& pos : _positions) {
      const auto lambda = _domain.x2Lambda(pos);
      const auto p = static_cast<size_t>(lambda[0] * m_partitionsCount);
      m_partition.push_back(p);

      const double px = (lambda[0] - p * partitionSize) / partitionSize;
      const auto x = static_cast<uint_fast32_t>(px * m_mortonRes);
      const auto y = static_cast<uint_fast32_t>(lambda[1] * m_mortonRes);
      const auto hash = libmorton::morton2D_64_encode(x, y);
      m_mortonHash.push_back(hash);
    }
  }

  void sortNodes(nodes& _nodes)
  {
    if (m_partition.size() != _nodes.size()) {
      throw std::runtime_error(
          "PartitionGenerator::sortGlobal: partition size does not match node "
          "size");
    }

    _nodes.reorder(
        ranges::view::zip(m_partition, m_mortonHash),
        [](const auto& _a, const auto& _b)
        {
          if (std::get<0>(_a) < std::get<0>(_b))  // Sort by partition index
            return true;
          if (std::get<0>(_a) > std::get<0>(_b))  // Sort by partition index
            return false;
          return std::get<1>(_a) < std::get<1>(_b);  // Sort by morton hash};
        });
  };

  void sortBonds(bonded::bonds& _bonds, const nodes& _nodes)
  {
    const NodeMap nodeMap = _nodes.getNodeMap();
    _bonds.remap(nodeMap);
    _bonds.flipSrcDst();

    checkPasses(_bonds);

    _bonds.reorder(_bonds.getBonds(),
                   [](const auto& _a, const auto& _b)
                   {
                     if (_a.src < _b.src)  // Sort by source node
                       return true;
                     if (_a.src > _b.src)  // Sort by source node
                       return false;
                     return _a.dst < _b.dst;  // Sort by destination node
                   });
  }

  auto getPasses() const -> size_t { return m_passes; }

  auto generatePartitions(const nodes& _nodes,
                          const bonded::bonds& _bonds) -> Partitions
  {
    Partitions partitions;
    partitions.reserve(m_partitionsCount);

    const auto partitionSize = _nodes.size() / m_partitionsCount;
    for (size_t i = 0; i < m_partitionsCount; ++i) {
      auto nodefirst = std::find(m_partition.begin(), m_partition.end(), i);
      auto nodelast =
          std::find(m_partition.rbegin(), m_partition.rend(), i).base();

      const size_t nodesStart = std::distance(m_partition.begin(), nodefirst);
      const size_t nodesEnd = std::distance(m_partition.begin(), nodelast);

      auto bondfirst = std::find_if(_bonds.getBonds().begin(),
                                    _bonds.getBonds().end(),
                                    [&](const auto& _bond)
                                    { return m_partition[_bond.src] == i; });
      auto bondlast = std::find_if(_bonds.getBonds().rbegin(),
                                   _bonds.getBonds().rend(),
                                   [&](const auto& _bond)
                                   { return m_partition[_bond.src] == i; })
                          .base();

      const size_t bondsStart =
          std::distance(_bonds.getBonds().begin(), bondfirst);
      const size_t bondsEnd =
          std::distance(_bonds.getBonds().begin(), bondlast);

      partitions.emplace_back(i, nodesStart, nodesEnd, bondsStart, bondsEnd);
    }

    return partitions;
  }

  void checkPasses(const bonded::bonds& _bonds)
  {
    for (const auto& bond : _bonds.getBonds()) {
      const size_t srcPar = m_partition[bond.src];
      const size_t dstPar = m_partition[bond.dst];

      const size_t srcPass = srcPar % m_passes;
      const size_t dstPass = dstPar % m_passes;

      // bonds need to be in the same partition or in seperate passes
      if ((srcPar == srcPar) || (srcPass != dstPass))
      {  // In the same partition
        continue;
      }
      throw std::runtime_error(
          "PartitionGenerator::sortBonds: bonds need to be "
          "in the same partition or in seperate passes");
    }
  }

private:
  std::size_t m_partitionsCount = 1;

  std::vector<size_t> m_partition;
  std::vector<std::uint_fast64_t> m_mortonHash;

  const size_t m_mortonRes = config::partition::mortonRes;
  size_t m_passes = 1;  //= config::partition::passes;
};

}  // namespace partition
}  // namespace networkV4