#include "Tensor2.hpp"

#include <cmath>
#include <cstddef>

template<class T>
tensor2<T>::tensor2()
    : xx(T())
    , xy(T())
    , yx(T())
    , yy(T())
{
}

template<class T>
tensor2<T>::tensor2(T _xx, T _xy, T _yx, T _yy)
    : xx(_xx)
    , xy(_xy)
    , yx(_yx)
    , yy(_yy)
{
}

template<class T>
tensor2<T>::tensor2(const tensor2& _t)
    : xx(_t.xx)
    , xy(_t.xy)
    , yx(_t.yx)
    , yy(_t.yy)
{
}

template<class T>
auto tensor2<T>::operator=(const tensor2& _t) -> tensor2&
{
    xx = _t.xx;
    xy = _t.xy;
    yx = _t.yx;
    yy = _t.yy;
    return *this;
}

template<class T>
auto tensor2<T>::operator+(const T& _v) const -> const tensor2
{
    return tensor2(xx + _v, xy + _v, yx + _v, yy + _v);
}

template<class T>
auto tensor2<T>::operator-(const T& _v) const -> const tensor2
{
    return tensor2(xx - _v, xy - _v, yx - _v, yy - _v);
}

template<class T>
auto tensor2<T>::operator*(const T& _s) const -> const tensor2
{
    return tensor2(xx * _s, xy * _s, yx * _s, yy * _s);
}

template<class T>
auto tensor2<T>::operator/(const T& _s) const -> const tensor2
{
    return tensor2(xx / _s, xy / _s, yx / _s, yy / _s);
}

template<class T>
auto tensor2<T>::operator+(const tensor2& _v) const -> const tensor2
{
    return tensor2(xx + _v.xx, xy + _v.xy, yx + _v.yx, yy + _v.yy);
}

template<class T>
auto tensor2<T>::operator-(const tensor2& _v) const -> const tensor2
{
    return tensor2(xx - _v.xx, xy - _v.xy, yx - _v.yx, yy - _v.yy);
}

template<class T>
auto tensor2<T>::operator*(const tensor2& _v) const -> const tensor2
{
    return tensor2(xx * _v.xx, xy * _v.xy, yx * _v.yx, yy * _v.yy);
}

template<class T>
auto tensor2<T>::operator/(const tensor2& _v) const -> const tensor2
{
    return tensor2(xx / _v.xx, xy / _v.xy, yx / _v.yx, yy / _v.yy);
}

template<class T>
auto tensor2<T>::operator+=(const T& _v) -> tensor2&
{
    xx += _v;
    xy += _v;
    yx += _v;
    yy += _v;
    return *this;
}

template<class T>
auto tensor2<T>::operator-=(const T& _v) -> tensor2&
{
    xx -= _v;
    xy -= _v;
    yx -= _v;
    yy -= _v;
    return *this;
}

template<class T>
auto tensor2<T>::operator*=(const T& _s) -> tensor2&
{
    xx *= _s;
    xy *= _s;
    yx *= _s;
    yy *= _s;
    return *this;
}

template<class T>
auto tensor2<T>::operator/=(const T& _s) -> tensor2&
{
    xx /= _s;
    xy /= _s;
    yx /= _s;
    yy /= _s;
    return *this;
}

template<class T>
auto tensor2<T>::operator+=(const tensor2& _v) -> tensor2&
{
    xx += _v.xx;
    xy += _v.xy;
    yx += _v.yx;
    yy += _v.yy;
    return *this;
}

template<class T>
auto tensor2<T>::operator-=(const tensor2& _v) -> tensor2&
{
    xx -= _v.xx;
    xy -= _v.xy;
    yx -= _v.yx;
    yy -= _v.yy;
    return *this;
}

template<class T>
auto tensor2<T>::operator*=(const tensor2& _v) -> tensor2&
{
    xx *= _v.xx;
    xy *= _v.xy;
    yx *= _v.yx;
    yy *= _v.yy;
    return *this;
}

template<class T>
auto tensor2<T>::operator/=(const tensor2& _v) -> tensor2&
{
    xx /= _v.xx;
    xy /= _v.xy;
    yx /= _v.yx;
    yy /= _v.yy;
    return *this;
}

template<class T>
void tensor2<T>::set(const T& _xx, const T& _xy, const T& _yx, const T& _yy)
{
    xx = _xx;
    xy = _xy;
    yx = _yx;
    yy = _yy;
}

template<class T>
void tensor2<T>::set(const tensor2& _t)
{
    xx = _t.xx;
    xy = _t.xy;
    yx = _t.yx;
    yy = _t.yy;
}

template<class T>
void tensor2<T>::set(const T& _v)
{
    xx = _v;
    xy = _v;
    yx = _v;
    yy = _v;
}

template<typename T>
tensor2<T> outer(const vec2<T>& _v1, const vec2<T>& _v2)
{
    return tensor2<T>(_v1.x * _v2.x, _v1.x * _v2.y, _v1.y * _v2.x, _v1.y * _v2.y);
}