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
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"
#include "OMP.hpp"

#if defined(_OPENMP)

networkV4::partition::Partitions threadPartitions;
size_t passes;

void networkV4::network::computeForces(bool _evalBreak)
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();

  for (size_t pass = 0; pass < passes; ++pass) {
    auto passParts =
        OMP::threadPartitions | ranges::views::drop(pass) | ranges::views::stride(OMP::passes);
    computePass(passParts, _evalBreak);
  }
}

void networkV4::network::computePass(auto _parts, bool _evalBreak)
{
#  pragma omp parallel for reduction(+ : m_energy) reduction(+ : m_stresses) \
      reduction(+ : m_breakQueue)
  for (const auto part : _parts) {
    for (auto&& [bond, type, brk] :
         ranges::views::zip(
             m_bonds.getBonds(), m_bonds.getTypes(), m_bonds.getBreaks())
             | ranges::views::slice(part.bondStart(), part.bondEnd()))
    {
      const auto& pos1 = m_nodes.positions()[bond.src];
      const auto& pos2 = m_nodes.positions()[bond.dst];
      const auto dist = m_box.minDist(pos1, pos2);

      if (_evalBreak) {
        const bool broken = bonded::visitBreak(brk, dist);
        if (broken) {
          m_breakQueue.emplace_back(bond.index, type, brk);

          type = Forces::VirtualBond {};
          brk = BreakTypes::None {};

          m_bondTags.addTag(bond.index, m_tags.get("broken"));
        }
      }

      const auto force = bonded::visitForce(type, dist);
      if (force) {
        m_nodes.forces()[bond.src] -= force.value();
        m_nodes.forces()[bond.dst] += force.value();

        const auto& tags = m_bondTags.getTagIds(bond.index);
        const auto stress = Utils::outer(force.value(), dist) / m_box.area();
        m_stresses.distribute(stress, tags);
      }

      const auto energy = bonded::visitEnergy(type, dist);
      if (energy) {
        m_energy += energy.value();
      }
    }
  }
}
#endif