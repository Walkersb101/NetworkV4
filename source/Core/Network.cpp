#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <ranges>
#include <vector>

#include "Network.hpp"

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Misc/Tensor2.hpp"
#include "Misc/Vec2.hpp"

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

auto networkV4::network::getBox() -> box&
{
  return m_box;
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

void networkV4::network::applyBond(const bonded::LocalBondInfo& _binfo,
                                   bool _evalBreak)
{
  const auto& type = m_bonds.getMap().getType(_binfo.index);
  const auto& breakType = m_bonds.getMap().getBreak(_binfo.index);
  const auto& tags = m_bonds.getMap().getTagIds(_binfo.index);

  const Utils::vec2d dist = m_box.minDist(m_nodes.positions()[_binfo.src],
                                          m_nodes.positions()[_binfo.dst]);

  if (_evalBreak && evalBreak(breakType, dist)) {
    m_breakQueue.push_back(std::make_pair(_binfo.index, breakType));  // TODO: needs to be reduced over openmp

    m_bonds.getMap().setType(_binfo.index, Forces::VirtualBond {});
    m_bonds.getMap().setBreak(_binfo.index, BreakTypes::None {});
    m_bonds.getMap().addTag(_binfo.index, m_tags.get("broken"));
    return;
  }

  const auto force = evalForce(type, dist);

  if (force) {
    m_nodes.forces()[_binfo.src] += force.value();
    m_nodes.forces()[_binfo.dst] -= force.value();

    m_stresses.distribute(Utils::outer(dist, force.value()) / m_box.area(),
                          tags);  // TODO: needs to be reduced over openmp
  }

  auto energy = evalEnergy(type, dist);
  if (energy) {
    m_energy += energy.value();  // TODO: needs to be reduced over openmp
  }
}

void networkV4::network::computeForces()
{
  m_energy = 0.0;
  m_stresses.zero();
  m_nodes.zeroForce();
  for (const auto& bond : m_bonds.getLocal()) {
    applyBond(bond);
  }
}

auto networkV4::network::computeEnergy() -> double
{
  m_energy = 0.0;
  for (const auto& bond : m_bonds.getLocal()) {
    const auto& type = m_bonds.getMap().getType(bond.index);
    const auto& dist = m_box.minDist(m_nodes.positions()[bond.src],
                                     m_nodes.positions()[bond.dst]);
    auto energy = evalEnergy(type, dist);
    if (energy) {
      m_energy += energy.value();
    }
  }
  return m_energy;
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