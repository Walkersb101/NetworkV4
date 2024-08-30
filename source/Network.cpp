#include <algorithm>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <unordered_map>
#include <vector>

#include "Network.hpp"

#include "Bonds.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"

networkV4::network::network()
    : m_nodes()
    , m_bonds()
    , m_stresses()
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

  _file.read(reinterpret_cast<char*>(&m_domain.x), sizeof(double));
  m_domain.y = m_domain.x;
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

  _file.read(reinterpret_cast<char*>(&m_domain), sizeof(vec2d));
  _file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));
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

auto networkV4::network::getStresses()
    -> std::unordered_map<bondType, tensor2d>&
{
  return m_stresses;
}

auto networkV4::network::getStresses() const
    -> const std::unordered_map<bondType, tensor2d>&
{
  return m_stresses;
}

auto networkV4::network::getShearStrain() const -> double
{
  return m_shearStrain;
}

auto networkV4::network::getElongationStrain() const -> double
{
  return m_elongationStrain;
}

auto networkV4::network::getShearStrain() -> double&
{
  return m_shearStrain;
}

auto networkV4::network::getElongationStrain() -> double&
{
  return m_elongationStrain;
}

auto networkV4::network::getRestSize() const -> vec2d
{
  return m_restSize;
}
auto networkV4::network::getDomain() const -> vec2d
{
  return m_domain;
}

auto networkV4::network::getRestSize() -> vec2d&
{
  return m_restSize;
}
auto networkV4::network::getDomain() -> vec2d&
{
  return m_domain;
}

void networkV4::network::applyBond(
    const bond& _bond,
    const nodes& _nodes,
    std::vector<vec2d>& _forces,
    std::unordered_map<bondType, tensor2d>& _stresses)
{
  if (!_bond.connected()) {
    return;
  }
  const vec2d dist =
      minDist(_nodes.position(_bond.src()), _nodes.position(_bond.dst()));
  const double r = dist.length();
  const vec2d force = dist * _bond.force(r) / r;
  _forces[_bond.src()] += force;
  _forces[_bond.dst()] -= force;

  _stresses[_bond.type()] += outer(force, dist);
}

void networkV4::network::applyBond(const bond& _bond)
{
  applyBond(_bond, m_nodes, m_nodes.forces(), m_stresses);
}

void networkV4::network::computeForces()
{
  clearStresses();
  getNodes().clearForces();
#pragma omp parallel for reduction(vec_vec2d_plus : m_nodes.forces()) \
    reduction(stresses_plus : m_stresses)
  for (const auto& _bond : m_bonds) {
    applyBond(_bond, m_nodes, m_nodes.forces(), m_stresses);
  }
  normalizeStresses();
}

auto networkV4::network::minDist(const vec2d& _pos1,
                                 const vec2d& _pos2) const -> vec2d
{
  const vec2d hDomain = m_domain * 0.5;
  const vec2d yShift = vec2d(m_domain.x * m_shearStrain, m_domain.y);
  vec2d dist = _pos2 - _pos1;
  while (dist.y > hDomain.y) {
    dist -= yShift;
  }
  while (dist.y <= -hDomain.y) {
    dist += yShift;
  }
  while (dist.x > hDomain.x) {
    dist.x -= m_domain.x;
  }
  while (dist.x < -hDomain.x) {
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
  const vec2d yShift = vec2d(m_domain.x * m_shearStrain, m_domain.y);
  while (_pos.y > m_domain.y) {
    _pos -= yShift;
  }
  while (_pos.y <= 0.0) {
    _pos += yShift;
  }
  while (_pos.x > m_domain.x) {
    _pos.x -= m_domain.x;
  }
  while (_pos.x < 0.0) {
    _pos.x += m_domain.x;
  }
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
}

void networkV4::network::elongate(double _step)
{
  m_elongationStrain += _step;
  const vec2d oldDomain = m_domain;
  m_domain = m_restSize
      * vec2d(1 / (1.0 + m_elongationStrain), 1.0 + m_elongationStrain);

  const vec2d scale = m_domain / oldDomain;
  std::transform(m_nodes.positions().begin(),
                 m_nodes.positions().end(),
                 m_nodes.positions().begin(),
                 [&](const vec2d& _pos) { return _pos * scale; });
}

void networkV4::network::initStresses()
{
  m_stresses.clear();
  for (const auto& bond : m_bonds) {
    if (m_stresses.find(bond.type()) == m_stresses.end()) {
      m_stresses[bond.type()] = tensor2d();
    }
  }
}

void networkV4::network::clearStresses()
{
  for (auto& stress : m_stresses) {
    stress.second = tensor2d();
  }
}

void networkV4::network::normalizeStresses()
{
  for (auto& stress : m_stresses) {
    stress.second /= m_domain.x * m_domain.y;
  }
}

auto networkV4::network::getGlobalStress() const -> tensor2d
{
  tensor2d stress;
  for (const auto& s : m_stresses) {
    stress += s.second;
  }
  return stress;
}

auto networkV4::network::bondStrain(const bond& _bond) const -> double
{
  const vec2d dist =
      minDist(m_nodes.position(_bond.src()), m_nodes.position(_bond.dst()));
  const double r = dist.length();
  return (r - _bond.naturalLength()) / _bond.naturalLength();
}

auto networkV4::network::bondEnergy(const bond& _bond) const -> double
{
  const double strain = bondStrain(_bond);
  return 0.5 * _bond.mu() * strain * strain;
}