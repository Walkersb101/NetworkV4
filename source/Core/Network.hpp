#pragma once

#include <deque>
#include <filesystem>
#include <vector>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Core/Stresses.hpp"
#include "Core/box.hpp"
#include "Misc/Enums.hpp"
#include "Misc/Tags/TagMap.hpp"
#include "Misc/Tags/TagStorage.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{

using bondQueue =
    std::deque<std::tuple<size_t, bonded::bondTypes, bonded::breakTypes>>;

class network
{
public:
  network() = delete;
  network(const box& _box, const size_t _N, const size_t _B);

public:
  auto getNodes()
      -> nodes&;  // TODO: make so Input and intergrators can access nodes
  auto getNodes() const -> const nodes&;

  auto getBonds() -> bonded::bonds&;  // TODO: make so Input and intergrators
                                      // can access nodes
  auto getBonds() const -> const bonded::bonds&;

  auto getEnergy() const -> const double;

  auto getRestBox() const -> const box;
  auto getBox() const -> const box;

  auto getTags() -> Utils::Tags::tagMap&;
  auto getTags() const -> const Utils::Tags::tagMap&;

  auto getBondTags() -> Utils::Tags::tagStorage&;
  auto getBondTags() const -> const Utils::Tags::tagStorage&;

  auto getNodeTags() -> Utils::Tags::tagStorage&;
  auto getNodeTags() const -> const Utils::Tags::tagStorage&;

public:
  auto getStresses()
      -> stresses&;  // TODO: make so Input and intergrators can access nodes
  auto getStresses() const -> const stresses&;

  auto getBreakQueue()
      -> bondQueue&;  // TODO: make so Input and intergrators can access nodes
  auto getBreakQueue() const -> const bondQueue&;

public:
  double getShearStrain() const;
  auto getElongationStrain() const -> Utils::vec2d;

  void shear(double _step);
  void setBox(const box& _box);

public:
  void computeForces();
  auto computeEnergy() -> double;

private:
  void applyBond(const bonded::BondInfo& _binfo,
                 bonded::bondTypes& _type,
                 bonded::breakTypes& _break,
                 bool _evalBreak = false);

private:
  box m_box;
  box m_restbox;

  double m_energy;
  stresses m_stresses;

  nodes m_nodes;
  bonded::bonds m_bonds;

  bondQueue m_breakQueue;

  Utils::Tags::tagMap m_tags;
  Utils::Tags::tagStorage m_bondTags;
  Utils::Tags::tagStorage m_nodeTags;
};

}  // namespace networkV4
