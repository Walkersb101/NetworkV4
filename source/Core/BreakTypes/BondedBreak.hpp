#pragma once

#include <variant>

#include "Core/BreakTypes/None.hpp"
#include "Core/BreakTypes/Strain.hpp"

namespace networkV4
{
namespace bonded
{

using breakTypes = std::variant<BreakTypes::None, BreakTypes::StrainBreak>;

inline auto visitBreak(const networkV4::bonded::breakTypes& _break,
                       const Utils::vec2d& _dist) -> bool
{
  return std::visit([_dist](const auto& _break) -> bool
                    { return _break.checkBreak(_dist); },
                    _break);
}

inline auto visitThreshold(const networkV4::bonded::breakTypes& _break,
                           const Utils::vec2d& _dist) -> std::optional<double>
{
  return std::visit([_dist](const auto& _break) -> std::optional<double>
                    { return _break.thresholdData(_dist); },
                    _break);
}

inline auto visitData(const networkV4::bonded::breakTypes& _break,
                           const Utils::vec2d& _dist) -> std::optional<double>
{
  return std::visit([_dist](const auto& _break) -> std::optional<double>
                    { return _break.data(_dist); },
                    _break);
}

}  // namespace bonded
}  // namespace networkV4