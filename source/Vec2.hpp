#pragma once

#include <cmath>
#include <algorithm>

template<class T>
struct vec2
{
  T x, y;
  vec2()
      : x(0)
      , y(0)
  {
  }
  vec2(T _x, T _y)
      : x(_x)
      , y(_y)
  {
  }

public:
  auto operator=(const vec2& _v) -> vec2&
  {
    x = _v.x;
    y = _v.y;
    return *this;
  }

  auto operator+(const T& _v) const -> const vec2
  {
    return vec2(x + _v, y + _v);
  }
  auto operator-(const T& _v) const -> const vec2
  {
    return vec2(x - _v, y - _v);
  }
  auto operator*(const T& _s) const -> const vec2
  {
    return vec2(x * _s, y * _s);
  }
  auto operator/(const T& _s) const -> const vec2
  {
    return vec2(x / _s, y / _s);
  }

  auto operator+(const vec2& _v) const -> const vec2
  {
    return vec2(x + _v.x, y + _v.y);
  }
  auto operator-(const vec2& _v) const -> const vec2
  {
    return vec2(x - _v.x, y - _v.y);
  }
  auto operator*(const vec2& _v) const -> const vec2
  {
    return vec2(x * _v.x, y * _v.y);
  }
  auto operator/(const vec2& _v) const -> const vec2
  {
    return vec2(x / _v.x, y / _v.y);
  }

  auto operator+=(const T& _v) -> vec2&
  {
    x += _v;
    y += _v;
    return *this;
  }
  auto operator-=(const T& _v) -> vec2&
  {
    x -= _v;
    y -= _v;
    return *this;
  }
  auto operator*=(const T& _s) -> vec2&
  {
    x *= _s;
    y *= _s;
    return *this;
  }
  auto operator/=(const T& _s) -> vec2&
  {
    x /= _s;
    y /= _s;
    return *this;
  }

  auto operator+=(const vec2& _v) -> vec2&
  {
    x += _v.x;
    y += _v.y;
    return *this;
  }
  auto operator-=(const vec2& _v) -> vec2&
  {
    x -= _v.x;
    y -= _v.y;
    return *this;
  }
  auto operator*=(const vec2& _v) -> vec2&
  {
    x *= _v.x;
    y *= _v.y;
    return *this;
  }
  auto operator/=(const vec2& _v) -> vec2&
  {
    x /= _v.x;
    y /= _v.y;
    return *this;
  }

  void set(const T& _x, const T& _y)
  {
    x = _x;
    y = _y;
  }
  void set(const vec2& _v)
  {
    x = _v.x;
    y = _v.y;
  }

  void normalize()
  {
    T l = length();
    x /= l;
    y /= l;
  }
  auto length() const -> T { return std::sqrt(x * x + y * y); }
  auto lengthSquared() const -> T { return x * x + y * y; }

  auto abs() const -> vec2 { return vec2(std::abs(x), std::abs(y)); }

  auto dot(const vec2& _v) const -> T { return x * _v.x + y * _v.y; }
  auto cross(const vec2& _v) const -> T { return x * _v.y - y * _v.x; }

  auto max() const -> T { return std::max(x, y); }
  auto min() const -> T { return std::min(x, y); }
  auto max(const vec2& _v) const -> vec2
  {
    return vec2(std::max(x, _v.x), std::max(y, _v.y));
  }

  auto min(const vec2& _v) const -> vec2
  {
    return vec2(std::min(x, _v.x), std::min(y, _v.y));
  }

  auto clamp(const vec2& _min, const vec2& _max) const -> vec2
  {
    return vec2(std::clamp(x, _min.x, _max.x), std::clamp(y, _min.y, _max.y));
  }
};

using vec2f = vec2<float>;
using vec2d = vec2<double>;