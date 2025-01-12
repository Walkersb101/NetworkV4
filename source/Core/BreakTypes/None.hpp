#pragma once

#include <optional>

#include "Misc/Math/Vector.hpp"

namespace networkV4
{
namespace BreakTypes
{
class None
{
public:
  None() {};

public:
  bool checkBreak(const Utils::Math::vec2d& _r) const { return false; }
  std::optional<double> thresholdData(const Utils::Math::vec2d& _r) const
  {
    return {};
  }
  std::optional<double> data(const Utils::Math::vec2d& _r) const { return {}; }
};
}  // namespace BreakTypes
}  // namespace networkV4