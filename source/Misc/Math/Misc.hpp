#pragma once

#include <cmath>

#include "Misc/Math/Matrix.hpp"
#include "Misc/Math/Vector.hpp"

namespace Utils
{
namespace Math
{

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> tensorProduct(const Vector<T, N>& _a, const Vector<T, M>& _b)
{
  Matrix<T, N, M> ret;
  auto ret_it = ret.begin();

  for (auto a : _a) {
    for (auto b : _b) {
      *(ret_it++) = a * b;
    }
  }
  return ret;
}

template<typename T, std::size_t N>
auto xdoty(const std::vector<Vector<T, N>>& _x, const std::vector<Vector<T, N>>& _y) -> T
{
  return std::inner_product(_x.begin(), _x.end(), _y.begin(), T(0));
}

}  // namespace Math
}  // namespace Utils