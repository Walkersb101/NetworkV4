#pragma once

#include <cstdint>
#include <ranges>
#include <vector>

#include "Misc/Config.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{
class nodes
{
public:
  nodes(bool _overdamped = false, bool _hasMass = true)
      : m_hasMass(_hasMass)
  {
  }

public:
  void clear();
  void reserve(size_t _size);

private:
  std::vector<size_t> m_globalIndices;

  std::vector<Utils::vec2d> m_positions;
  std::vector<Utils::vec2d> m_forces;
  std::vector<Utils::vec2d> m_velocities;

  const bool m_hasMass;
  std::vector<double> m_masses;
};
}  // namespace networkV4