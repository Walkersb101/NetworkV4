#pragma once

#include <optional>

#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace BreakTypes
{
class None
{
public:
  None() {};

public:
  bool checkBreak(const Utils::vec2d& _r) const {
  return false;
}
  std::optional<double> thresholdData(const Utils::vec2d& _r) const {
  return {};
}
};
}  // namespace BreakTypes
}  // namespace networkV4