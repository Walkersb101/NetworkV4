#pragma once

#include <algorithm>
#include <array>
#include <utility>
#include <cmath>
#include <numeric>

namespace Utils
{
namespace Math
{

template<typename T, std::size_t N>
class Vector : public std::array<T, N>
{
  using Base = std::array<T, N>;

public:
  using Base::Base;
  Vector() = default;
  Vector(Vector const&) = default;
  Vector& operator=(Vector const&) = default;

public:
  using Base::at;
  using Base::operator[];
  using Base::front;
  using Base::back;
  using Base::data;
  using Base::begin;
  using Base::cbegin;
  using Base::end;
  using Base::cend;
  using Base::empty;
  using Base::size;
  using Base::max_size;
  using Base::fill;
  using array_type = std::array<T, N>;

public:
  void swap(Vector& rhs) { std::swap_ranges(begin(), end(), rhs.begin()); }

public:
  T norm2() const { return (*this) * (*this); }
  T norm() const { return std::sqrt(norm2()); }

  Vector& normalize()
  {
    auto const l = norm();
    if (l > T(0)) {
      std::transform(this->begin(),
                     this->end(),
                     this->begin(),
                     [l](T& val) { return val / l; });
    }
    return *this;
  }

  Vector normalized() const { return (*this) / (*this).norm(); }

  Vector abs()
  {
    Vector<T, N> ret;
    std::transform(this->begin(), this->end(), ret.begin(), [](T const& val) {
      return std::abs(val);
    });
    return ret;
  }

  T max() { return *std::max_element(this->begin(), this->end()); }
  T min() { return *std::min_element(this->begin(), this->end()); }
};

template<class T>
using Vector3 = Vector<T, 3>;

template<std::size_t N>
using VectorXd = Vector<double, N>;
using Vector2d = VectorXd<2>;
using Vector3d = VectorXd<3>;
using Vector4d = VectorXd<4>;
using Vector6d = VectorXd<6>;
using Vector9d = VectorXd<9>;

template<std::size_t N>
using VectorXf = Vector<float, N>;
using Vector3f = VectorXf<3>;

template<std::size_t N>
using VectorXi = Vector<int, N>;
using Vector3i = VectorXi<3>;

// ---- Helper functions ----

template<std::size_t N, typename T, typename U, typename Op>
auto binary_op(Vector<T, N> const& a, Vector<U, N> const& b, Op op)
{
  using R = decltype(op(std::declval<T>(), std::declval<U>()));
  Vector<R, N> ret;

  std::transform(
      std::begin(a), std::end(a), std::begin(b), std::begin(ret), op);

  return ret;
}

template <std::size_t N, typename T, typename Op>
Vector<T, N> &binary_op_assign(Vector<T, N> &a, Vector<T, N> const &b, Op op) {
  std::transform(std::begin(a), std::end(a), std::begin(b), std::begin(a), op);
  return a;
}

template<std::size_t N, typename T, typename Op>
constexpr bool all_of(Vector<T, N> const& a, Vector<T, N> const& b, Op op)
{
    return std::all_of(std::begin(a), std::end(a), std::begin(b), [&op](const T& ai, const T& bi) {
    return static_cast<bool>(op(ai, bi)); });
}

template<class T>
struct is_vector : std::false_type
{
};
template<class T, std::size_t N>
struct is_vector<Vector<T, N>> : std::true_type
{
};

// ---- Comparison operators ----

template<std::size_t N, typename T>
constexpr bool operator<(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return all_of(a, b, std::less<T>());
}

template<std::size_t N, typename T>
constexpr bool operator>(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return all_of(a, b, std::greater<T>());
}

template<std::size_t N, typename T>
constexpr bool operator<=(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return all_of(a, b, std::less_equal<T>());
}

template<std::size_t N, typename T>
constexpr bool operator>=(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return all_of(a, b, std::greater_equal<T>());
}

template<std::size_t N, typename T>
constexpr bool operator==(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return all_of(a, b, std::equal_to<T>());
}

template<std::size_t N, typename T>
constexpr bool operator!=(Vector<T, N> const& a, Vector<T, N> const& b)
{
  return not(a == b);
}

// ---- addition and subtraction ----

template<std::size_t N, typename T, typename U>
auto operator+(Vector<T, N> const& a, Vector<U, N> const& b)
{
  return binary_op(a, b, std::plus<>());
}

template<std::size_t N, typename T>
Vector<T, N>& operator+=(Vector<T, N>& a, Vector<T, N> const& b)
{
  return binary_op_assign(a, b, std::plus<T>());
}

template<std::size_t N, typename T, typename U>
auto operator-(Vector<T, N> const& a, Vector<U, N> const& b)
{
  return binary_op(a, b, std::minus<>());
}

template<std::size_t N, typename T>
Vector<T, N> operator-(Vector<T, N> const& a)
{
  Vector<T, N> ret;

  std::transform(std::begin(a), std::end(a), std::begin(ret), std::negate<T>());

  return ret;
}

template<std::size_t N, typename T>
Vector<T, N>& operator-=(Vector<T, N>& a, Vector<T, N> const& b)
{
  return binary_op_assign(a, b, std::minus<T>());
}

// ---- scalar multiplication and division ----

template<std::size_t N,
         typename T,
         class U,
         std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(U const& a, Vector<T, N> const& b)
{
  using R = decltype(a * std::declval<T>());
  Vector<R, N> ret;

  std::transform(std::begin(b),
                 std::end(b),
                 std::begin(ret),
                 [a](T const& val) { return a * val; });

  return ret;
}

template<std::size_t N,
         typename T,
         class U,
         std::enable_if_t<std::is_arithmetic_v<U>, bool> = true>
auto operator*(Vector<T, N> const& b, U const& a)
{
  using R = decltype(std::declval<T>() * a);
  Vector<R, N> ret;

  std::transform(std::begin(b),
                 std::end(b),
                 std::begin(ret),
                 [a](T const& val) { return a * val; });

  return ret;
}

template<std::size_t N, typename T>
Vector<T, N>& operator*=(Vector<T, N>& b, T const& a)
{
  std::transform(std::begin(b),
                 std::end(b),
                 std::begin(b),
                 [a](T const& val) { return a * val; });
  return b;
}

template<std::size_t N, typename T>
Vector<T, N> operator/(Vector<T, N> const& a, T const& b)
{
  Vector<T, N> ret;

  std::transform(std::begin(a),
                 std::end(a),
                 ret.begin(),
                 [b](T const& val) { return val / b; });
  return ret;
}

template<std::size_t N, typename T>
Vector<T, N> operator/(T const& a, Vector<T, N> const& b)
{
  Vector<T, N> ret;

  std::transform(std::begin(b),
                 std::end(b),
                 ret.begin(),
                 [a](T const& val) { return a / val; });
  return ret;
}

template<std::size_t N, typename T>
Vector<T, N>& operator/=(Vector<T, N>& a, T const& b)
{
  std::transform(std::begin(a),
                 std::end(a),
                 std::begin(a),
                 [b](T const& val) { return val / b; });
  return a;
}

// ---- dot product ----

template<
    std::size_t N,
    typename T,
    class U,
    class = std::enable_if_t<not(is_vector<T>::value or is_vector<U>::value)>>
auto operator*(Vector<T, N> const& a, Vector<U, N> const& b)
{
  using std::declval;
  using R = decltype(std::declval<T>() * std::declval<U>());

  return std::inner_product(std::begin(a), std::end(a), std::begin(b), R {});
}

}  // namespace math
}  // namespace Utils
