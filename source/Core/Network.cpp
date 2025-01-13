#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <vector>

#include "Network.hpp"

#include <range/v3/algorithm.hpp>
#include <range/v3/view/zip.hpp>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Misc/Math/Vector.hpp"
#include "Misc/Tensor2.hpp"

networkV4::network::network(const box& _box, const size_t _N, const size_t _B)
    : m_box(_box)
    , m_restbox(_box)
    , m_energy(0.0)
    , m_stresses()
    , m_nodes(_N)
    , m_bonds(_B)
    , m_breakQueue()
    , m_tags()
{
  // add default tags
  m_tags.add("broken");
}

auto networkV4::network::getNodes() -> nodes&
{
  return m_nodes;
}

auto networkV4::network::getNodes() const -> const nodes&
{
  return m_nodes;
}

auto networkV4::network::getBonds() -> bonded::bonds&
{
  return m_bonds;
}

auto networkV4::network::getBonds() const -> const bonded::bonds&
{
  return m_bonds;
}

auto networkV4::network::getEnergy() const -> const double
{
  return m_energy;
}

auto networkV4::network::getRestBox() const -> const box
{
  return m_restbox;
}

auto networkV4::network::getBox() const -> const box
{
  return m_box;
}

auto networkV4::network::getTags() -> Utils::Tags::tagMap&
{
  return m_tags;
}

auto networkV4::network::getTags() const -> const Utils::Tags::tagMap&
{
  return m_tags;
}

auto networkV4::network::getBondTags() -> Utils::Tags::tagStorage&
{
  return m_bondTags;
}

auto networkV4::network::getBondTags() const -> const Utils::Tags::tagStorage&
{
  return m_bondTags;
}

auto networkV4::network::getNodeTags() -> Utils::Tags::tagStorage&
{
  return m_nodeTags;
}

auto networkV4::network::getNodeTags() const -> const Utils::Tags::tagStorage&
{
  return m_nodeTags;
}

auto networkV4::network::getStresses() -> stresses&
{
  return m_stresses;
}

auto networkV4::network::getStresses() const -> const stresses&
{
  return m_stresses;
}

auto networkV4::network::getBreakQueue() -> bondQueue&
{
  return m_breakQueue;
}

auto networkV4::network::getBreakQueue() const -> const bondQueue&
{
  return m_breakQueue;
}

double networkV4::network::getShearStrain() const
{
  return m_box.shearStrain();
}

auto networkV4::network::getElongationStrain() const -> Utils::Math::vec2d
{
  return broadcast(m_box.getDomain() - m_restbox.getDomain(),
                   m_restbox.getDomain(),
                   std::divides<>());
}

void networkV4::network::shear(double _step)
{
  double dxy = _step * m_box.getLy();
  m_box.setxy(m_box.getxy() + dxy);
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](const Utils::Math::vec2d& _pos) {
                   return _pos + Utils::Math::vec2d({_step * _pos[1], 0.0});
                 });
}

void networkV4::network::setBox(const box& _box)
{
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](Utils::Math::vec2d& _pos)
                 {
                   const auto plambda = m_box.x2Lambda(_pos);
                   return _box.lambda2x(plambda);
                 });
  m_box = _box;
}

void networkV4::network::wrapNodes()
{
  for (auto& pos : m_nodes.positions()) {
    pos = m_box.wrapPosition(pos);
  }
}

#if not defined(_OPENMP)
template <bool _evalBreak, bool _evalStress>
void networkV4::network::computeForces()
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();

  const auto& positions = m_nodes.positions();

  for (auto&& [bond, type, brk, tags] : ranges::views::zip(m_bonds.getBonds(),
                                                           m_bonds.getTypes(),
                                                           m_bonds.getBreaks(),
                                                           m_bonds.getTags()))
  {
    const auto& pos1 = positions[bond.src];
    const auto& pos2 = positions[bond.dst];
    const auto dist = m_box.minDist(pos1, pos2);

    if constexpr (_evalBreak) {
      evalBreak(dist, bond, type, brk, tags);
    }

    const auto force = bonded::visitForce(type, dist);
    if (force) {
      applyforce(bond, dist, force.value(), tags, _evalStress);
    }

    const auto energy = bonded::visitEnergy(type, dist);
    if (energy) {
      m_energy += energy.value();
    }
  }
}

// Explicit template instantiation
template void networkV4::network::computeForces<false, false>();
template void networkV4::network::computeForces<true, false>();
template void networkV4::network::computeForces<false, true>();
template void networkV4::network::computeForces<true, true>();
#endif

auto networkV4::network::computeEnergy() -> double
{
  m_energy = 0.0;
  for (const auto [bond, type] :
       ranges::views::zip(m_bonds.getBonds(), m_bonds.getTypes()))
  {
    const auto& pos1 = m_nodes.positions()[bond.src];
    const auto& pos2 = m_nodes.positions()[bond.dst];
    const auto dist = m_box.minDist(pos1, pos2);

    const auto energy = bonded::visitEnergy(type, dist);
    if (energy) {
      m_energy += energy.value();
    }
  }
  return m_energy;
}

void networkV4::network::computeBreaks()
{
  for (const auto& [bond, type, brk, tags] :
       ranges::views::zip(m_bonds.getBonds(),
                          m_bonds.getTypes(),
                          m_bonds.getBreaks(),
                          m_bonds.getTags()))
  {
    const auto& pos1 = m_nodes.positions()[bond.src];
    const auto& pos2 = m_nodes.positions()[bond.dst];
    const auto dist = m_box.minDist(pos1, pos2);

    evalBreak(dist, bond, type, brk, tags);
  }
}

void networkV4::network::evalBreak(const Utils::Math::vec2d& _dist,
                                   const bonded::BondInfo& _binfo,
                                   bonded::bondTypes& _type,
                                   bonded::breakTypes& _break,
                                   Utils::Tags::tagFlags& _tags)
{
  const bool broken = bonded::visitBreak(_break, _dist);
  if (broken) {
    m_breakQueue.emplace_back(_binfo, _type, _break, _tags);

    _type = Forces::VirtualBond {};
    _break = BreakTypes::None {};

    _tags.set(BROKEN_TAG_INDEX);
  }
}

void networkV4::network::applyforce(const bonded::BondInfo& _binfo,
                                    const Utils::Math::vec2d& _dist,
                                    const Utils::Math::vec2d& _force,
                                    const Utils::Tags::tagFlags& _tags,
                                    bool _evalStress)
{
  auto& forces = m_nodes.forces();
  forces[_binfo.src] -= _force;
  forces[_binfo.dst] += _force;

  if (_evalStress) {
    const auto stress =
        Utils::Math::tensorProduct(_force, _dist) * m_box.invArea();
    m_stresses.distribute(stress, _tags);
  }
}