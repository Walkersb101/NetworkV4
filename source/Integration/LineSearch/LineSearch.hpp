#pragma once

#include <cstdint>

namespace networkV4
{
namespace lineSearch
{
enum class lineSearchState : uint8_t
{
  success,
  DirectionNotDescent,
  zeroforce,
  zeroquad,
  zeroAlpha,
};

}  // namespace lineSearch
}  // namespace networkV4