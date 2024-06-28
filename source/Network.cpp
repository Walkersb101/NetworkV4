#include "Network.hpp"

#include <cmath>
#include <unordered_map>

#include "Bonds.hpp"

network::network::network()
    : m_nodes()
    , m_bonds()
{
}

network::network::~network() {}

auto network::network::getNodes() -> nodes&
{
  return m_nodes;
}

auto network::network::getNodes() const -> const nodes&
{
  return m_nodes;
}

auto network::network::getStresses() -> std::unordered_map<bondType, tensor2d>&
{
  return m_stresses;
}

void network::network::computeForces()
{
  clearStresses();
  for (const auto& b : m_bonds) {
    if (!b.connected()) {
      continue;
    }
    const vec2d dist =
        minDist(m_nodes.position(b.src()), m_nodes.position(b.dst()));
    const double r = dist.length();
    const vec2d force =  dist * b.force(r) / r;
    m_nodes.force(b.src()) += force;
    m_nodes.force(b.dst()) -= force;

    const tensor2d stress = outer(force, dist);
    m_stresses[b.type()] += stress;
  }
}

auto network::network::minDist(const vec2d& _pos1,
                               const vec2d& _pos2) const -> vec2d
{
  vec2d dist = _pos2 - _pos1;
  dist.x -= std::round(dist.y / m_domain.y) * m_domain.x * m_shearStrain;
  dist.x = std::remainder(dist.x, m_domain.x);
  dist.y = std::remainder(dist.y, m_domain.y);
  return dist;
  /*
    vec2d dist = _pos2 - _pos1;
    while (dist.y > m_domain.y / 2.0) {
      dist -= m_domain * vec2d(-m_shearStrain, 1.0);
    }

    while (dist.y < -m_domain.y / 2.0) {
      dist += m_domain * vec2d(m_shearStrain, 1.0);
    }

    while (dist.x > m_domain.x / 2.0) {
      dist.x -= m_domain.x;
    }

    while (dist.x < -m_domain.x / 2.0) {
      dist.x += m_domain.x;
    }
    return dist;
    */
}

void network::network::clearStresses()
{
  for (auto& stress : m_stresses) {
    stress.second = tensor2d();
  }
}

void network::network::normalizeStresses()
{
  for (auto& stress : m_stresses) {
    stress.second /= m_domain.x * m_domain.y;
  }
}
