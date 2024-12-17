#pragma once

#include <cmath>

#include "Vec2.hpp"

namespace Utils
{
template<class T>
struct tensor2
{
public:
  T xx, xy, yx, yy;
  tensor2()
      : xx(0)
      , xy(0)
      , yx(0)
      , yy(0)
  {
  }
  tensor2(T _xx, T _xy, T _yx, T _yy)
      : xx(_xx)
      , xy(_xy)
      , yx(_yx)
      , yy(_yy)
  {
  }
  tensor2(const tensor2& _t)
      : xx(_t.xx)
      , xy(_t.xy)
      , yx(_t.yx)
      , yy(_t.yy)
  {
  }

public:
  auto operator=(const tensor2& _t) -> tensor2&
  {
    xx = _t.xx;
    xy = _t.xy;
    yx = _t.yx;
    yy = _t.yy;
    return *this;
  }

  auto operator+(const T& _v) const -> const tensor2
  {
    return tensor2(xx + _v, xy + _v, yx + _v, yy + _v);
  }
  auto operator-(const T& _v) const -> const tensor2
  {
    return tensor2(xx - _v, xy - _v, yx - _v, yy - _v);
  }
  auto operator*(const T& _s) const -> const tensor2
  {
    return tensor2(xx * _s, xy * _s, yx * _s, yy * _s);
  }
  auto operator/(const T& _s) const -> const tensor2
  {
    return tensor2(xx / _s, xy / _s, yx / _s, yy / _s);
  }

  auto operator+(const tensor2& _v) const -> const tensor2
  {
    return tensor2(xx + _v.xx, xy + _v.xy, yx + _v.yx, yy + _v.yy);
  }
  auto operator-(const tensor2& _v) const -> const tensor2
  {
    return tensor2(xx - _v.xx, xy - _v.xy, yx - _v.yx, yy - _v.yy);
  }
  auto operator*(const tensor2& _v) const -> const tensor2
  {
    return tensor2(xx * _v.xx, xy * _v.xy, yx * _v.yx, yy * _v.yy);
  }
  auto operator/(const tensor2& _v) const -> const tensor2
  {
    return tensor2(xx / _v.xx, xy / _v.xy, yx / _v.yx, yy / _v.yy);
  }

  auto operator+=(const T& _v) -> tensor2&
  {
    xx += _v;
    xy += _v;
    yx += _v;
    yy += _v;
    return *this;
  }
  auto operator-=(const T& _v) -> tensor2&
  {
    xx -= _v;
    xy -= _v;
    yx -= _v;
    yy -= _v;
    return *this;
  }
  auto operator*=(const T& _s) -> tensor2&
  {
    xx *= _s;
    xy *= _s;
    yx *= _s;
    yy *= _s;
    return *this;
  }
  auto operator/=(const T& _s) -> tensor2&
  {
    xx /= _s;
    xy /= _s;
    yx /= _s;
    yy /= _s;
    return *this;
  }

  auto operator+=(const tensor2& _v) -> tensor2&
  {
    xx += _v.xx;
    xy += _v.xy;
    yx += _v.yx;
    yy += _v.yy;
    return *this;
  }
  auto operator-=(const tensor2& _v) -> tensor2&
  {
    xx -= _v.xx;
    xy -= _v.xy;
    yx -= _v.yx;
    yy -= _v.yy;
    return *this;
  }
  auto operator*=(const tensor2& _v) -> tensor2&
  {
    xx *= _v.xx;
    xy *= _v.xy;
    yx *= _v.yx;
    yy *= _v.yy;
    return *this;
  }
  auto operator/=(const tensor2& _v) -> tensor2&
  {
    xx /= _v.xx;
    xy /= _v.xy;
    yx /= _v.yx;
    yy /= _v.yy;
    return *this;
  }

  void set(const T& _xx, const T& _xy, const T& _yx, const T& _yy)
  {
    xx = _xx;
    xy = _xy;
    yx = _yx;
    yy = _yy;
  }
  void set(const tensor2& _t)
  {
    xx = _t.xx;
    xy = _t.xy;
    yx = _t.yx;
    yy = _t.yy;
  }
  void set(const T& _v)
  {
    xx = _v;
    xy = _v;
    yx = _v;
    yy = _v;
  }
};

template<typename T>
auto outer(const vec2<T>& _v1, const vec2<T>& _v2) -> tensor2<T> {
  return tensor2<T>(_v1.x * _v2.x, _v1.x * _v2.y, _v1.y * _v2.x, _v1.y * _v2.y);
}

using tensor2f = tensor2<float>;
using tensor2d = tensor2<double>;
} // namespace Utils