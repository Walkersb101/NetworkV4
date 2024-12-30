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
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

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

auto networkV4::network::getElongationStrain() const -> Utils::vec2d
{
  return (m_box.getDomain() - m_restbox.getDomain()) / m_restbox.getDomain();
}

void networkV4::network::shear(double _step)
{
  double dxy = _step * m_box.getLy();
  m_box.setxy(m_box.getxy() + dxy);
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](const Utils::vec2d& _pos)
                 { return _pos + Utils::vec2d(_step * _pos.y, 0.0); });
}

void networkV4::network::setBox(const box& _box)
{
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](Utils::vec2d& _pos)
                 {
                   const auto plambda = m_box.x2Lambda(_pos);
                   return _box.lambda2x(plambda);
                 });
  m_box = _box;
}

inline auto evalForce(const networkV4::bonded::bondTypes& _bond,
                      const Utils::vec2d& _dist) -> std::optional<Utils::vec2d>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<Utils::vec2d>
                    { return _bond.force(_dist); },
                    _bond);
}

inline auto evalEnergy(const networkV4::bonded::bondTypes& _bond,
                       const Utils::vec2d& _dist) -> std::optional<double>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<double>
                    { return _bond.energy(_dist); },
                    _bond);
}

inline auto evalBreak(const networkV4::bonded::breakTypes& _break,
                      const Utils::vec2d& _dist) -> bool
{
  return std::visit([_dist](const auto& _break) -> bool
                    { return _break.checkBreak(_dist); },
                    _break);
}

void networkV4::network::computeForces()
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();
  for (const auto& [bond, type, brk] :
       ranges::views::zip(m_bonds.getBonds(),
                          m_bonds.getTypes(),
                          m_bonds.getBreaks()))
    {
    //applyBond(bond);
  }
}

auto networkV4::network::computeEnergy() -> double
{
  m_energy = 0.0;
  for (const auto [bond, type] :
       ranges::views::zip(m_bonds.getBonds(), m_bonds.getTypes()))
  {
    const auto& pos1 = m_nodes.positions()[bond.src];
    const auto& pos2 = m_nodes.positions()[bond.dst];
    const auto dist = m_box.minDist(pos1, pos2);

    m_energy += evalEnergy(type, dist).value_or(0.0);
  }
  return m_energy;
}

void networkV4::network::applyBond(const bonded::BondInfo& _binfo,
                                   bonded::bondTypes& _type,
                                   bonded::breakTypes& _break,
                                   bool _evalBreak)
{
  const auto& pos1 = m_nodes.positions()[_binfo.src];
  const auto& pos2 = m_nodes.positions()[_binfo.dst];
  const auto dist = m_box.minDist(pos1, pos2);

  if (_evalBreak && evalBreak(_break, dist)) {
    m_breakQueue.emplace_back(_binfo.index,
                              _type,
                              _break);  // TODO: needs to be reduced over openmp

    _type = Forces::VirtualBond {};
    _break = BreakTypes::None {};

    m_bondTags.addTag(_binfo.index, m_tags.get("broken"));
    return;
  }

  const auto force = evalForce(_type, dist);

  if (force) {
    m_nodes.forces()[_binfo.src] += force.value();
    m_nodes.forces()[_binfo.dst] -= force.value();

    const auto& tags = m_bondTags.getTagIds(_binfo.index);
    const auto stress = Utils::outer(dist, force.value()) / m_box.area();
    m_stresses.distribute(stress,
                          tags);  // TODO: needs to be reduced over openmp
  }

  m_energy += evalEnergy(_type, dist)
                  .value_or(0.0);  // TODO: needs to be reduced over openmp
}

/*
void networkV4::network::shear(double _step)
{
  m_shearStrain += _step;
  const double offset = _step * m_domain.y * 0.5;
  std::transform(
      m_nodes.positions().begin(),
      m_nodes.positions().end(),
      m_nodes.positions().begin(),
      [&](const vec2d& _pos)
      { return wrapedPosition(_pos + vec2d((_step * _pos.y) - offset, 0.0));
}); m_yshift = vec2d(m_domain.x * m_shearStrain, m_domain.y);
}

void networkV4::network::elongate(double _step)
{
  m_elongationStrain += _step;
  const vec2d oldDomain = m_domain;
  setDomain(m_restSize
            * vec2d(1 / (1.0 + m_elongationStrain), 1.0 +
m_elongationStrain));

  const vec2d scale = m_domain / oldDomain;
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](const vec2d& _pos) { return _pos * scale; });
}
*/