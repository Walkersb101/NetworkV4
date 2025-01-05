#pragma once

#include <span>
#include <vector>

#include "Misc/Vec2.hpp"

namespace Utils
{

template<typename T>
auto spanView(const std::vector<vec2<T>>& _vec) -> std::span<const T>
{
  return std::span<const T>(reinterpret_cast<const T*>(_vec.data()),
                            _vec.size() * 2);
}

template<typename T>
auto spanView(std::vector<vec2<T>>& _vec) -> std::span<T>
{
  return std::span<T>(reinterpret_cast<T*>(_vec.data()), _vec.size() * 2);
}

}  // namespace Utils