#include <algorithm>

#include "Core/Nodes.hpp"

networkV4::nodes::nodes(bool _overdamped, bool _hasMass)
    : m_overdamped(_overdamped)
    , m_hasMass(_hasMass)
{
}

void networkV4::nodes::clear()
{
  m_positions.clear();
  m_forces.clear();
  m_velocities.clear();
  m_masses.clear();
}

void networkV4::nodes::reserve(size_t _size)
{
  m_positions.reserve(_size);
  m_forces.reserve(_size);

  if (!m_overdamped) {
    m_velocities.reserve(_size);
  }
  if (m_hasMass) {
    m_masses.reserve(_size);
  }
}