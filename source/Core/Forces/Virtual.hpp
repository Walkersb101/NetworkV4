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
  VirtualBond();

public:
  std::optional<Utils::vec2d> force(const Utils::vec2d& _r) const;
  std::optional<double> energy(const Utils::vec2d& _r) const;
};

VirtualBond::VirtualBond() {}

inline std::optional<Utils::vec2d> VirtualBond::force(
    const Utils::vec2d& _dx) const
{
  return {};
}

inline std::optional<double> VirtualBond::energy(const Utils::vec2d& _dx) const
{
  return {};
}
}  // namespace Forces
}  // namespace networkV4