#pragma once

#include <span>
#include <vector>

#include "Misc/Math/Vector.hpp"

namespace Utils
{

template<typename T, std::size_t N>
auto spanView(const std::vector<Math::Vector<T, N>>& _vec) -> std::span<const T>
{
  return std::span<const T>(reinterpret_cast<const T*>(_vec.data()), _vec.size() * N);
}

template<typename T, std::size_t N>
auto spanView(std::vector<Math::Vector<T, N>>& _vec) -> std::span<T>
{
  return std::span<T>(reinterpret_cast<T*>(_vec.data()), _vec.size() * N);
}

template<typename T>
inline auto sign(const T& _val) -> T
{
  return _val < 0 ? -1 : 1;
}

inline auto rms(const std::vector<Math::vec2d>& _vec) -> double
{
  double sum = 0.0;
  for (const auto& val : _vec) {
    sum += val.norm2();
  }
  return std::sqrt(sum / _vec.size());
}

inline auto maxComp(const std::vector<Math::vec2d>& _vec) -> double
{
  double max = 0.0;
  for (const auto& val : _vec) {
    max = std::max(max, val.abs().max());
  }
  return max;
}

}  // namespace Utils