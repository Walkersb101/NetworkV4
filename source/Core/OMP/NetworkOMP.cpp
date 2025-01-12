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
#include "Misc/Math/Vector.hpp"
#include "OMP.hpp"

#if defined(_OPENMP)

networkV4::partition::Partitions networkV4::OMP::threadPartitions;
size_t networkV4::OMP::passes;

void networkV4::network::computeForces(bool _evalBreak)
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();

  for (size_t pass = 0; pass < OMP::passes; ++pass) {
    auto passParts = OMP::threadPartitions | ranges::views::drop(pass)
        | ranges::views::stride(OMP::passes);
    computePass(passParts, _evalBreak);
  }
}

void networkV4::network::computePass(auto _parts, bool _evalBreak)
{
  const auto& bonds = m_bonds.getBonds();
  auto& types = m_bonds.getTypes();
  auto& breaks = m_bonds.getBreaks();
  auto& tags = m_bonds.getTags();

  const auto& positions = m_nodes.positions();
  auto& forces = m_nodes.forces();

#  pragma omp parallel for reduction(+ : m_energy) reduction(+ : m_stresses) \
      reduction(+ : m_breakQueue) \
      num_threads(OMP::threadPartitions.size() / 2) schedule(static, 1)
  for (const auto part : _parts) {
    for (size_t i = part.bondStart(); i < part.bondEnd(); i++) {
      const auto& bond = bonds[i];
      auto& type = types[i];
      auto& brk = breaks[i];
      auto& bTags = tags[i];

      const auto& pos1 = positions[bond.src];
      const auto& pos2 = positions[bond.dst];
      const auto dist = m_box.minDist(pos1, pos2);

      if (_evalBreak) {
        const bool broken = bonded::visitBreak(brk, dist);
        if (broken) {
          m_breakQueue.emplace_back(bond, type, brk, bTags);

          type = Forces::VirtualBond {};
          brk = BreakTypes::None {};

          bTags.set(BROKEN_TAG_INDEX);
        }
      }

      const auto force = bonded::visitForce(type, dist);
      if (force) {
        const auto& f = force.value();
        forces[bond.src] -= f;
        forces[bond.dst] += f;

        const auto stress = Utils::Math::tensorProduct(f, dist) * m_box.invArea();
        m_stresses.distribute(stress, bTags);
      }

      const auto energy = bonded::visitEnergy(type, dist);
      if (energy) {
        m_energy += energy.value();
      }
    }
  }
}
#endif