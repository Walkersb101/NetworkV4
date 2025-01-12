#pragma once

#include <variant>

#include "Core/Forces/Harmonic.hpp"
#include "Core/Forces/Virtual.hpp"

namespace networkV4
{
namespace bonded
{

using bondTypes = std::variant<Forces::VirtualBond, Forces::HarmonicBond>;

inline auto bondName(const bondTypes& _bond) -> std::string
{
  return std::visit(
      [](auto&& arg) -> std::string
      {
        using T = std::decay_t<decltype(arg)>;
        if constexpr (std::is_same_v<T, Forces::VirtualBond>) {
          return "VirtualBond";
        } else if constexpr (std::is_same_v<T, Forces::HarmonicBond>) {
          return "HarmonicBond";
        } else {
          throw("bondName: unknown bond type");
        }
      },
      _bond);
}

inline auto visitForce(const networkV4::bonded::bondTypes& _bond,
                       const Utils::Math::vec2d& _dist) -> std::optional<Utils::Math::vec2d>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<Utils::Math::vec2d>
                    { return _bond.force(_dist); },
                    _bond);
}

inline auto visitEnergy(const networkV4::bonded::bondTypes& _bond,
                        const Utils::Math::vec2d& _dist) -> std::optional<double>
{
  return std::visit([_dist](const auto& _bond) -> std::optional<double>
                    { return _bond.energy(_dist); },
                    _bond);
}

}  // namespace bonded
}  // namespace networkV4