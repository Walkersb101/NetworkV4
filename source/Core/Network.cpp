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
#include "Core/Tensor2.hpp"
#include "Core/Vec2.hpp"
#include "Misc/EnumMap.hpp"

networkV4::network::network()
    : m_nodes()
    , m_bonds()
    , m_stresses()
    , m_energy(0.0)
    , m_domain()
    , m_restSize()
    , m_shearStrain(0.0)
    , m_elongationStrain(0.0)
{
}

networkV4::network::~network() {}

void networkV4::network::loadFromBinV1(std::ifstream& _file, double _lambda)
{
  std::size_t N;
  std::size_t B;
  _file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  _file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  double domainX;
  _file.read(reinterpret_cast<char*>(&domainX), sizeof(double));
  setDomain(vec2d(domainX, domainX));
  _file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));
  m_restSize = m_domain;

  m_nodes.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    vec2d pos;
    _file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
    m_nodes.addNode(pos);
  }

  m_bonds.resize(B);
  for (auto& b : m_bonds) {
    bool connected;
    std::size_t index1;
    std::size_t index2;
    double naturalLength;
    double constant;

    _file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    _file.read(reinterpret_cast<char*>(&constant), sizeof(double));

    b = bond(index1,
             index2,
             connected,
             bondType::single,
             naturalLength,
             constant,
             _lambda);
  }
  initStresses();
}

void networkV4::network::loadFromBinV2(std::ifstream& _file)
{
  std::size_t N;
  std::size_t B;
  _file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  _file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  vec2d domain;  
  _file.read(reinterpret_cast<char*>(&domain), sizeof(vec2d));
  _file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));
  setDomain(domain);
  m_restSize = m_domain;

  m_nodes.reserve(N);
  for (size_t i = 0; i < N; ++i) {
    vec2d pos;
    _file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
    m_nodes.addNode(pos);
  }

  m_bonds.reserve(B);
  for (size_t i = 0; i < B; ++i) {
    std::size_t index1;
    std::size_t index2;
    bool connected;
    bool matrix;
    double naturalLength;
    double constant;
    double lambda;

    _file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    _file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&matrix), sizeof(bool));
    _file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    _file.read(reinterpret_cast<char*>(&constant), sizeof(double));
    _file.read(reinterpret_cast<char*>(&lambda), sizeof(double));

    if (index1 >= N || index2 >= N) {
      throw std::runtime_error("Invalid bond indices: " + std::to_string(index1)
                               + ", " + std::to_string(index2));
    }

    const bondType type = matrix ? bondType::matrix : bondType::sacrificial;
    m_bonds.add(
        bond(index1, index2, connected, type, naturalLength, constant, lambda));
  }
  initStresses();
}

auto networkV4::network::getNodes() -> nodes&
{
  return m_nodes;
}

auto networkV4::network::getNodes() const -> const nodes&
{
  return m_nodes;
}

auto networkV4::network::getBonds() -> bonds&
{
  return m_bonds;
}

auto networkV4::network::getBonds() const -> const bonds&
{
  return m_bonds;
}

auto networkV4::network::getStresses() -> stressMap&
{
  return m_stresses;
}

auto networkV4::network::getStresses() const -> const stressMap&
{
  return m_stresses;
}

auto networkV4::network::getEnergy() const -> const double
{
  return m_energy;
}

auto networkV4::network::getShearStrain() const -> const double
{
  return m_shearStrain;
}

auto networkV4::network::getElongationStrain() const -> const double
{
  return m_elongationStrain;
}

auto networkV4::network::getRestSize() const -> const vec2d
{
  return m_restSize;
}
auto networkV4::network::getDomain() const -> const vec2d
{
  return m_domain;
}

void networkV4::network::setDomain(const vec2d& _domain)
{
  m_domain = _domain;
  m_halfDomain = m_domain * 0.5;
  m_yshift = vec2d(m_domain.x * m_shearStrain, m_domain.y);
}

void networkV4::network::applyBond(const bond& _bond,
                                   const nodes& _nodes,
                                   std::vector<vec2d>& _forces,
                                   stressMap& _stresses,
                                   double& _energy)
{
  // if (!_bond.connected()) {
  //   return;
  // }
  const vec2d dist =
      minDist(_nodes.position(_bond.src()), _nodes.position(_bond.dst()));
  const double length = dist.length();
  const vec2d force = dist * (_bond.force(length) / length);
  _forces[_bond.src()] += force;
  _forces[_bond.dst()] -= force;

  _energy += _bond.energy(length);

  _stresses[_bond.type()] += outer(force, dist);
}

void networkV4::network::computeForces()
{
  m_energy = 0.0;
  clearStresses();
  getNodes().clearForces();
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = m_bonds | std::views::filter(connected);
  for (const auto& bond : connectedList) {
    applyBond(bond, m_nodes, m_nodes.forces(), m_stresses, m_energy);
  }
  normalizeStresses();
}

auto networkV4::network::computeEnergy() -> double
{
  m_energy = 0.0;
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = m_bonds | std::views::filter(connected);
  for (const auto& bond : connectedList) {
    const vec2d dist =
        minDist(m_nodes.position(bond.src()), m_nodes.position(bond.dst()));
    m_energy += bond.energy(dist.length());
  }
  return m_energy;
}

auto networkV4::network::minDist(const vec2d& _pos1,
                                 const vec2d& _pos2) const -> vec2d
{
  vec2d dist = _pos2 - _pos1;
  if (m_shearStrain != 0.0) {
    while (std::abs(dist.y) > m_halfDomain.y) {
      if (dist.y > 0.0)
        dist -= m_yshift;
      else
        dist += m_yshift;
    }
  } else {
    while (std::abs(dist.y) > m_halfDomain.y) {
      if (dist.y > 0.0)
        dist.y -= m_domain.y;
      else
        dist.y += m_domain.y;
    }
  }
  while (std::abs(dist.x) > m_halfDomain.x) {
    if (dist.x > 0.0)
      dist.x -= m_domain.x;
    else
      dist.x += m_domain.x;
  }
  // vec2d dist2 = _pos2 - _pos1;
  // dist2.x -= std::round(dist2.y / m_domain.y) * m_domain.x * m_shearStrain;
  // distBC(dist2.x, m_domain.x);
  // distBC(dist2.y, m_domain.y);
  return dist;
}

void networkV4::network::wrapPosition(vec2d& _pos) const
{
  if (m_shearStrain != 0.0) {
    while (std::abs(_pos.y) > m_domain.y) {
      if (_pos.y > 0.0)
        _pos -= m_yshift;
      else
        _pos += m_yshift;
    }
  } else {
    while (std::abs(_pos.y) > m_domain.y) {
      if (_pos.y > 0.0)
        _pos.y -= m_domain.y;
      else
        _pos.y += m_domain.y;
    }
  }
  while (std::abs(_pos.x) > m_domain.x) {
    if (_pos.x > 0.0)
      _pos.x -= m_domain.x;
    else
      _pos.x += m_domain.x;
  }

  //while (_pos.y > m_domain.y)
  //  _pos -= m_yshift;
  //while (_pos.y <= 0.0)
  //  _pos += m_yshift;
  //while (_pos.x > m_domain.x)
  //  _pos.x -= m_domain.x;
  //while (_pos.x < 0.0)
  //  _pos.x += m_domain.x;

  //_pos.x -= std::floor(_pos.y / m_domain.y) * m_domain.x * m_shearStrain;
  // BC(_pos.x, m_domain.x);
  // BC(_pos.y, m_domain.y);
}

auto networkV4::network::wrapedPosition(const vec2d& _pos) const -> vec2d
{
  vec2d pos = _pos;
  wrapPosition(pos);
  return pos;
}

void networkV4::network::shear(double _step)
{
  m_shearStrain += _step;
  const double offset = _step * m_domain.y * 0.5;
  std::transform(
      m_nodes.positions().begin(),
      m_nodes.positions().end(),
      m_nodes.positions().begin(),
      [&](const vec2d& _pos)
      { return wrapedPosition(_pos + vec2d((_step * _pos.y) - offset, 0.0)); });
  m_yshift = vec2d(m_domain.x * m_shearStrain, m_domain.y);
}

void networkV4::network::elongate(double _step)
{
  m_elongationStrain += _step;
  const vec2d oldDomain = m_domain;
  setDomain(m_restSize
            * vec2d(1 / (1.0 + m_elongationStrain), 1.0 + m_elongationStrain));

  const vec2d scale = m_domain / oldDomain;
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](const vec2d& _pos) { return _pos * scale; });
}

void networkV4::network::initStresses()
{
  m_stresses = stressMap();
}

void networkV4::network::clearStresses()
{
  m_stresses = stressMap();
}

void networkV4::network::normalizeStresses()
{
  for (auto& stress : m_stresses) {
    stress /= m_domain.x * m_domain.y;
  }
}

auto networkV4::network::getGlobalStress() const -> tensor2d
{
  return std::accumulate(m_stresses.begin(),
                         m_stresses.end(),
                         tensor2d(),
                         [](const tensor2d& _sum, const auto& _stress)
                         { return _sum + _stress; });
}

auto networkV4::network::bondStrain(const bond& _bond) const -> double
{
  const vec2d dist =
      minDist(m_nodes.position(_bond.src()), m_nodes.position(_bond.dst()));
  return (dist.length() - _bond.naturalLength()) / _bond.naturalLength();
}

auto networkV4::network::bondEnergy(const bond& _bond) const -> double
{
  const vec2d dist =
      minDist(m_nodes.position(_bond.src()), m_nodes.position(_bond.dst()));
  return 0.5 * _bond.mu() * std::pow(dist.length() - _bond.naturalLength(), 2) / _bond.naturalLength();
}