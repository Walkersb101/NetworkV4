#pragma once

#include <cmath>

#include "Misc/Math/Matrix.hpp"
#include "Misc/Math/Vector.hpp"

namespace Utils
{
namespace Math
{

template<typename T, std::size_t N, std::size_t M>
Matrix<T, N, M> tensorProduct(Vector<T, N> const& _a, Vector<T, M> const& _b)
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

}  // namespace Math
}  // namespace Utils