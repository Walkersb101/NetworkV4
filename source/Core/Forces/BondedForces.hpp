#pragma once

#include <variant>

#include "Core/Forces/Harmonic.hpp"
#include "Core/Forces/Virtual.hpp"

namespace networkV4
{
namespace bonded
{

using bondTypes = std::variant<Forces::VirtualBond, Forces::HarmonicBond>;

inline auto visitForce(const networkV4::bonded::bondTypes& _bond,
                      const Utils::vec2d& _dist) -> std::optional<Utils::vec2d>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<Utils::vec2d>
                    { return _bond.force(_dist); },
                    _bond);
}

inline auto visitEnergy(const networkV4::bonded::bondTypes& _bond,
                       const Utils::vec2d& _dist) -> std::optional<double>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<double>
                    { return _bond.energy(_dist); },
                    _bond);
}

}  // namespace bonded
}  // namespace networkV4