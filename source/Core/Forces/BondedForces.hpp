#pragma once

#include <variant>

#include "Core/Forces/Harmonic.hpp"
#include "Core/Forces/Virtual.hpp"

namespace networkV4
{
namespace bonded
{

using bondTypes = std::variant<Forces::VirtualBond, Forces::HarmonicBond>;

}  // namespace bonded
}  // namespace networkV4