#pragma once

#include <optional>

#include "Misc/Config.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace Forces
{
class VirtualBond
{
public:
  VirtualBond() {}

public:
  std::optional<Utils::vec2d> force(const Utils::vec2d& _dx) const
  {
    return {};
  }
  std::optional<double> energy(const Utils::vec2d& _dx) const { return {}; }
};
}  // namespace Forces
}  // namespace networkV4