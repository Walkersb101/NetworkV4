#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iostream>
#include <numeric>
#include <vector>

#include "Core/Network.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view/slice.hpp>
#include <range/v3/view/zip.hpp>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Misc/Math/Tensor2.hpp"
#include "Misc/Math/Vector.hpp"
#include "OMP.hpp"

#if defined(_OPENMP)

networkV4::partition::Partitions networkV4::OMP::threadPartitions;
size_t networkV4::OMP::passes;
std::vector<networkV4::stresses> networkV4::OMP::localStresses;
std::vector<networkV4::bondQueue> networkV4::OMP::localBreaks;

template<bool _evalBreak, bool _evalStress>
void networkV4::network::computeForces()
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();

  for (size_t pass = 0; pass < OMP::passes; ++pass) {
    auto passParts = OMP::threadPartitions | ranges::views::drop(pass)
        | ranges::views::stride(OMP::passes);
    computePass<_evalBreak, _evalStress>(passParts);
  }
}

template<bool _evalBreak, bool _evalStress>
void networkV4::network::computePass(auto _parts)
{
  const auto& bonds = m_bonds.getBonds();
  auto& types = m_bonds.getTypes();
  auto& breaks = m_bonds.getBreaks();
  auto& tags = m_bonds.getTags();

  const auto& positions = m_nodes.positions();
  auto& forces = m_nodes.forces();

#  pragma omp parallel for reduction(+ : m_energy) num_threads(_parts.size()) \
      schedule(static, 1)
  for (const auto part : _parts) {
    size_t threadID = omp_get_thread_num();
    auto& localStresses = OMP::localStresses[threadID];
    auto& localBreaks = OMP::localBreaks[threadID];

    if constexpr (_evalStress) {
      localStresses.zero();
      localStresses.init(m_stresses.getInitilised());
    }

    if constexpr (_evalBreak) {
      localBreaks.clear();
    }

    for (size_t i = part.bondStart(); i < part.bondEnd(); i++) {
      const auto& bond = bonds[i];
      auto& type = types[i];
      auto& brk = breaks[i];
      auto& bTags = tags[i];

      if constexpr (!_evalBreak
                    && std::is_same_v<Forces::VirtualBond, decltype(type)>)
      {
        continue;
      }

      const auto& pos1 = positions[bond.src];
      const auto& pos2 = positions[bond.dst];
      const auto dist = m_box.minDist(pos1, pos2);

      if constexpr (_evalBreak) {
        const bool broken = bonded::visitBreak(brk, dist);
        if (broken) {
          localBreaks.emplace_back(bond, type, brk, bTags);

          type = Forces::VirtualBond {};
          brk = BreakTypes::None {};

          bTags.set(BROKEN_TAG_INDEX);
        }
      }

      const auto force = bonded::visitForce(type, dist);
      if (force) {
        const auto& f = force.value();
        forces[bond.src] += f;
        forces[bond.dst] -= f;

        if constexpr (_evalStress) {
          const auto stress =
              Utils::Math::tensorProduct(f, -dist) * m_box.invArea();
          localStresses.distribute(stress, bTags);
        }
      }

      const auto energy = bonded::visitEnergy(type, dist);
      if (energy) {
        m_energy += energy.value();
      }
    }

    if constexpr (_evalStress) {
#  pragma omp critical
      merge(m_stresses, localStresses);
    }

    if constexpr (_evalBreak) {
#  pragma omp critical
      merge(m_breakQueue, localBreaks);
    }
  }
}

// Explicit template instantiation
template void networkV4::network::computeForces<false, false>();
template void networkV4::network::computeForces<true, false>();
template void networkV4::network::computeForces<false, true>();
template void networkV4::network::computeForces<true, true>();

template void networkV4::network::computePass<false, false>(
    decltype(OMP::threadPartitions));
template void networkV4::network::computePass<true, false>(
    decltype(OMP::threadPartitions));
template void networkV4::network::computePass<false, true>(
    decltype(OMP::threadPartitions));
template void networkV4::network::computePass<true, true>(
    decltype(OMP::threadPartitions));
#endif