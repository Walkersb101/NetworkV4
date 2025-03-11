#pragma once

#include <algorithm>
#include <array>
#include <cmath>
#include <numeric>
#include <utility>

#include "Misc/Math/Vector.hpp"

namespace Utils
{
namespace Math
{

template<typename T, std::size_t Rows, std::size_t Cols>
struct Matrix
{
  using container = Utils::Math::Vector<T, Cols * Rows>;
  using pointer = typename container::pointer;
  using const_pointer = typename container::const_pointer;
  using iterator = typename container::iterator;
  using const_iterator = typename container::const_iterator;
  using value_type = typename container::value_type;
  using reference = typename container::reference;
  using const_reference = typename container::const_reference;

  container m_data;

  Matrix() { std::fill(begin(), end(), T()); }
  Matrix(std::initializer_list<T> init_list)
  {
    assert(init_list.size() == Rows * Cols);
    std::copy(init_list.begin(), init_list.end(), begin());
  }
  Matrix(std::initializer_list<std::initializer_list<T>> init_list)
  {
    assert(init_list.size() == Rows);
    auto row_it = init_list.begin();
    for (std::size_t i = 0; i < Rows; ++i, ++row_it) {
      assert(row_it->size() == Cols);
      std::copy(row_it->begin(), row_it->end(), begin() + i * Cols);
    }
  }

  constexpr value_type operator()(std::size_t row, std::size_t col) const
  {
    assert(row < Rows);
    assert(col < Cols);
    return m_data[Cols * row + col];
  }

  constexpr pointer data() { return m_data.data(); }
  constexpr const_pointer data() const noexcept { return m_data.data(); }

  constexpr iterator begin() noexcept { return m_data.begin(); }
  constexpr const_iterator begin() const noexcept { return m_data.begin(); }
  constexpr iterator end() noexcept { return m_data.end(); }
  constexpr const_iterator end() const noexcept { return m_data.end(); }
};

template<typename T>
using mat22 = Matrix<T, 2, 2>;
using mat22d = mat22<double>;
using mat22f = mat22<float>;

template<typename T>
using mat33 = Matrix<T, 3, 3>;
using mat33d = mat33<double>;
using mat33f = mat33<float>;

template<typename T, std::size_t M, std::size_t N>
Utils::Math::Vector<T, M * N> flatten(Matrix<T, M, N> const& m)
{
  return Utils::Math::Vector<T, M * N>(m.begin(), m.end());
}

// ---- Helper functions ----

template<std::size_t Rows,
         std::size_t Cols,
         typename T,
         typename U,
         typename Op>
auto broadcast(Matrix<T, Rows, Cols> const& a,
               Matrix<U, Rows, Cols> const& b,
               Op op)
{
  using R = decltype(op(std::declval<T>(), std::declval<U>()));
  Matrix<R, Rows, Cols> ret;

  std::transform(a.begin(), a.end(), b.begin(), ret.begin(), op);

  return ret;
}

template<std::size_t Rows, std::size_t Cols, typename T, typename Op>
Matrix<T, Rows, Cols>& broadcast_assign(Matrix<T, Rows, Cols>& a,
                                        Matrix<T, Rows, Cols> const& b,
                                        Op op)
{
  std::transform(a.begin(), a.end(), b.begin(), a.begin(), op);
  return a;
}

template<std::size_t Rows, std::size_t Cols, typename T, typename Op>
constexpr bool all_of(Matrix<T, Rows, Cols> const& a,
                      Matrix<T, Rows, Cols> const& b,
                      Op op)
{
  return std::all_of(a.begin(),
                     a.end(),
                     b.begin(),
                     [&op](const T& ai, const T& bi)
                     { return static_cast<bool>(op(ai, bi)); });
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator<(Matrix<T, Rows, Cols> const& a,
                         Matrix<T, Rows, Cols> const& b)
{
  return all_of(a, b, std::less<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator>(Matrix<T, Rows, Cols> const& a,
                         Matrix<T, Rows, Cols> const& b)
{
  return all_of(a, b, std::greater<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator<=(Matrix<T, Rows, Cols> const& a,
                          Matrix<T, Rows, Cols> const& b)
{
  return all_of(a, b, std::less_equal<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator>=(Matrix<T, Rows, Cols> const& a,
                          Matrix<T, Rows, Cols> const& b)
{
  return all_of(a, b, std::greater_equal<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator==(Matrix<T, Rows, Cols> const& a,
                          Matrix<T, Rows, Cols> const& b)
{
  return all_of(a, b, std::equal_to<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
constexpr bool operator!=(Matrix<T, Rows, Cols> const& a,
                          Matrix<T, Rows, Cols> const& b)
{
  return not(a == b);
}

// ---- addition and subtraction ----

template<std::size_t Rows, std::size_t Cols, typename T, typename U>
auto operator+(Matrix<T, Rows, Cols> const& a, Matrix<U, Rows, Cols> const& b)
{
  return broadcast(a, b, std::plus<>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols>& operator+=(Matrix<T, Rows, Cols>& a,
                                  Matrix<T, Rows, Cols> const& b)
{
  return broadcast_assign(a, b, std::plus<T>());
}

template<std::size_t Rows, std::size_t Cols, typename T, typename U>
auto operator-(Matrix<T, Rows, Cols> const& a, Matrix<U, Rows, Cols> const& b)
{
  return broadcast(a, b, std::minus<>());
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols> operator-(Matrix<T, Rows, Cols> const& a)
{
  Matrix<T, Rows, Cols> ret;
  std::transform(a.begin(), a.end(), ret.begin(), std::negate<T>());
  return ret;
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols>& operator-=(Matrix<T, Rows, Cols>& a,
                                  Matrix<T, Rows, Cols> const& b)
{
  return broadcast_assign(a, b, std::minus<T>());
}

// Scalar multiplication and division
template<std::size_t Rows,
         std::size_t Cols,
         typename T,
         class U,
         std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(U const& a, Matrix<T, Rows, Cols> const& b)
{
  using R = decltype(a * std::declval<T>());
  Matrix<R, Rows, Cols> ret;
  std::transform(
      b.begin(), b.end(), ret.begin(), [a](T const& val) { return a * val; });
  return ret;
}

template<std::size_t Rows,
         std::size_t Cols,
         typename T,
         class U,
         std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(Matrix<T, Rows, Cols> const& b, U const& a)
{
  return a * b;
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols>& operator*=(Matrix<T, Rows, Cols>& b, T const& a)
{
  std::transform(
      b.begin(), b.end(), b.begin(), [a](T const& val) { return a * val; });
  return b;
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols> operator/(Matrix<T, Rows, Cols> const& a, T const& b)
{
  Matrix<T, Rows, Cols> ret;
  std::transform(
      a.begin(), a.end(), ret.begin(), [b](T const& val) { return val / b; });
  return ret;
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols> operator/(T const& a, Matrix<T, Rows, Cols> const& b)
{
  Matrix<T, Rows, Cols> ret;
  std::transform(
      b.begin(), b.end(), ret.begin(), [a](T const& val) { return a / val; });
  return ret;
}

template<std::size_t Rows, std::size_t Cols, typename T>
Matrix<T, Rows, Cols>& operator/=(Matrix<T, Rows, Cols>& a, T const& b)
{
  std::transform(
      a.begin(), a.end(), a.begin(), [b](T const& val) { return val / b; });
  return a;
}

}  // namespace Math
}  // namespace Utils