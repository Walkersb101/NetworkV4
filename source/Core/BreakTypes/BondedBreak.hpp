#pragma once

#include <variant>

#include "Core/BreakTypes/None.hpp"
#include "Core/BreakTypes/Strain.hpp"

namespace networkV4
{
namespace bonded
{
    
using breakTypes = std::variant<BreakTypes::None, BreakTypes::StrainBreak>;

}  // namespace bonded
}  // namespace networkV4