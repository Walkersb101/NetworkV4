#pragma once

#include <optional>

#include "Misc/Config.hpp"
#include "Misc/Math/Vector.hpp"

namespace networkV4
{
namespace Forces
{
class VirtualBond
{
public:
  VirtualBond() {}

public:
  std::optional<Utils::Math::vec2d> force(const Utils::Math::vec2d& _dx) const
  {
    return {};
  }
  std::optional<double> energy(const Utils::Math::vec2d& _dx) const { return {}; }
};
}  // namespace Forces
}  // namespace networkV4