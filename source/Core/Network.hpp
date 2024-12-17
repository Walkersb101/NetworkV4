#pragma once

#include <deque>
#include <filesystem>
#include <vector>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Core/box.hpp"
#include "Core/Stresses.hpp"
#include "Misc/Enums.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{
class network
{
public:
    network() = delete;
    network(const box& _box, const nodes& _nodes, const bonded::bonds& _bonds);

public:
  auto getNodes() -> nodes&;
  auto getNodes() const -> const nodes&;

  auto getBonds() -> bonded::bonds&;
  auto getBonds() const -> const bonded::bonds&;

  auto getEnergy() const -> const double;

  auto getRestBox() const -> const box;
  auto getBox() const -> const box;
  auto getBox() -> box&;

public:
  void computeForces();
  auto computeEnergy() -> double;

private:
  void applyBond(const bonded::LocalBondInfo& _binfo, bool _evalBreak = false);

private:
  box m_box;
  box m_restbox;

  double m_energy;
  stresses m_stresses;

  nodes m_nodes;
  bonded::bonds m_bonds;

  std::deque<std::pair<size_t, bonded::breakTypes>> m_breakQueue;

  tagMap m_tags;  // TODO: fill with defined tags (broken)
};

}  // namespace networkV4
