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
{
}

networkV4::network::~network() {}

void checkCanOpen(const std::filesystem::path& _filename)
{
  if (!std::filesystem::exists(_filename)) {
    throw std::runtime_error("File does not exist: " + _filename.string());
  }
  std::ifstream file(_filename, std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _filename.string());
  }
}

void networkV4::network::loadFromBinV1(const std::filesystem::path& _filename,
                                       double _lambda)
{
  checkCanOpen(_filename);
  std::ifstream file(_filename, std::ios::in | std::ios::binary);

  std::size_t N;
  std::size_t B;
  file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  file.read(reinterpret_cast<char*>(&m_domain.x), sizeof(double));
  m_domain.y = m_domain.x;
  file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));

  m_nodes.resize(N);
  for (vec2d& pos : m_nodes.positions()) {
    file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
  }

  m_bonds.resize(B);
  for (auto& b : m_bonds) {
    bool connected;
    std::size_t index1;
    std::size_t index2;
    double naturalLength;
    double constant;

    file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    file.read(reinterpret_cast<char*>(&constant), sizeof(double));

    b = bond(index1,
             index2,
             connected,
             bondType::single,
             naturalLength,
             constant,
             _lambda);
  }
}

void networkV4::network::loadFromBinV2(const std::filesystem::path& _filename)
{
  checkCanOpen(_filename);
  std::ifstream file(_filename, std::ios::in | std::ios::binary);

  std::size_t N;
  std::size_t B;
  file.read(reinterpret_cast<char*>(&N), sizeof(std::size_t));
  file.read(reinterpret_cast<char*>(&B), sizeof(std::size_t));

  file.read(reinterpret_cast<char*>(&m_domain.x), sizeof(double));
  file.read(reinterpret_cast<char*>(&m_domain.y), sizeof(double));
  file.read(reinterpret_cast<char*>(&m_shearStrain), sizeof(double));

  m_nodes.resize(N);
  for (vec2d& pos : m_nodes.positions()) {
    file.read(reinterpret_cast<char*>(&pos), sizeof(vec2d));
  }

  m_bonds.reserve(B);
  for (auto& b : m_bonds) {
    std::size_t index1;
    std::size_t index2;
    bool connected;
    bool matrix;
    double naturalLength;
    double constant;
    double lambda;

    file.read(reinterpret_cast<char*>(&index1), sizeof(std::size_t));
    file.read(reinterpret_cast<char*>(&index2), sizeof(std::size_t));
    file.read(reinterpret_cast<char*>(&connected), sizeof(bool));
    file.read(reinterpret_cast<char*>(&matrix), sizeof(bool));
    file.read(reinterpret_cast<char*>(&naturalLength), sizeof(double));
    file.read(reinterpret_cast<char*>(&constant), sizeof(double));
    file.read(reinterpret_cast<char*>(&lambda), sizeof(double));

    if (index1 >= N || index2 >= N) {
      throw std::runtime_error("Invalid bond indices: " + std::to_string(index1)
                               + ", " + std::to_string(index2));
    }

    const bondType type = matrix ? bondType::matrix : bondType::sacrificial;
    b = bond(index1, index2, connected, type, naturalLength, constant, lambda);
  }
}

void networkV4::network::loadFromBin(const std::filesystem::path& _filename,
                                     loadVersion _version,
                                     double _lambda)
{
  switch (_version) {
    case loadVersion::binV1:
      loadFromBinV1(_filename, _lambda);
      break;
    case loadVersion::binV2:
      loadFromBinV2(_filename);
      break;
    default:
      throw std::runtime_error(
          "Unknown version: "
          + std::to_string(static_cast<std::uint8_t>(_version)));
  }
}

auto networkV4::network::getNodes() -> nodes&
{
  return m_nodes;
}

auto networkV4::network::getNodes() const -> const nodes&
{
  return m_nodes;
}

auto networkV4::network::getStresses()
    -> std::unordered_map<bondType, tensor2d>&
{
  return m_stresses;
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
  vec2d dist = _pos2 - _pos1;
  dist.x -= std::round(dist.y / m_domain.y) * m_domain.x * m_shearStrain;
  dist.x = std::remainder(dist.x, m_domain.x);
  dist.y = std::remainder(dist.y, m_domain.y);
  return dist;
}

void networkV4::network::wrapPosition(vec2d& _pos) const
{
  _pos.x -= std::floor(_pos.y / m_domain.y) * m_domain.x * m_shearStrain;
  _pos.x = std::remainder(_pos.x, m_domain.x);
  _pos.y = std::remainder(_pos.y, m_domain.y);
}

auto networkV4::network::wrapedPosition(const vec2d& _pos) const -> vec2d
{
  vec2d pos = _pos;
  wrapPosition(pos);
  return pos;
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