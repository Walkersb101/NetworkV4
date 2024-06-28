#include "Vec2.hpp"

template<class T>
vec2<T>::vec2()
    : x(0)
    , y(0)
{
}

template<class T>
vec2<T>::vec2(T _x, T _y)
    : x(_x)
    , y(_y)
{
}

template<class T>
vec2<T>::vec2(const vec2& _v)
    : x(_v.x)
    , y(_v.y)
{
}

template<class T>
auto vec2<T>::operator=(const vec2& _v) -> vec2&
{
    x = _v.x;
    y = _v.y;
    return *this;
}

template<class T>
auto vec2<T>::operator+(const T& _v) const -> const vec2
{
    return vec2(x + _v, y + _v);
}

template<class T>
auto vec2<T>::operator-(const T& _v) const -> const vec2
{
    return vec2(x - _v, y - _v);
}

template<class T>
auto vec2<T>::operator*(const T& _s) const -> const vec2
{
    return vec2(x * _s, y * _s);
}

template<class T>
auto vec2<T>::operator/(const T& _s) const -> const vec2
{
    return vec2(x / _s, y / _s);
}

template<class T>
auto vec2<T>::operator+(const vec2& _v) const -> const vec2
{
    return vec2(x + _v.x, y + _v.y);
}

template<class T>
auto vec2<T>::operator-(const vec2& _v) const -> const vec2
{
    return vec2(x - _v.x, y - _v.y);
}

template<class T>
auto vec2<T>::operator*(const vec2& _v) const -> const vec2
{
    return vec2(x * _v.x, y * _v.y);
}

template<class T>
auto vec2<T>::operator/(const vec2& _v) const -> const vec2
{
    return vec2(x / _v.x, y / _v.y);
}

template<class T>
auto vec2<T>::operator+=(const T& _v) -> vec2&
{
    x += _v;
    y += _v;
    return *this;
}
template<class T>
auto vec2<T>::operator-=(const T& _v) -> vec2&
{
    x -= _v;
    y -= _v;
    return *this;
}

template<class T>
auto vec2<T>::operator*=(const T& _s) -> vec2&
{
    x *= _s;
    y *= _s;
    return *this;
}

template<class T>
auto vec2<T>::operator/=(const T& _s) -> vec2&
{
    x /= _s;
    y /= _s;
    return *this;
}

template<class T>
auto vec2<T>::operator+=(const vec2& _v) -> vec2&
{
    x += _v.x;
    y += _v.y;
    return *this;
}

template<class T>
auto vec2<T>::operator-=(const vec2& _v) -> vec2&
{
    x -= _v.x;
    y -= _v.y;
    return *this;
}

template<class T>
auto vec2<T>::operator*=(const vec2& _v) -> vec2&
{
    x *= _v.x;
    y *= _v.y;
    return *this;
}

template<class T>
auto vec2<T>::operator/=(const vec2& _v) -> vec2&
{
    x /= _v.x;
    y /= _v.y;
    return *this;
}

template<class T>
void vec2<T>::set(const T& _x, const T& _y)
{
    x = _x;
    y = _y;
}

template<class T>
void vec2<T>::set(const vec2& _v)
{
    x = _v.x;
    y = _v.y;
}

template<class T>
void vec2<T>::normalize()
{
    T l = length();
    x /= l;
    y /= l;
}

template<class T>
auto vec2<T>::length() const -> T
{
    return std::sqrt(x * x + y * y);
}

template<class T>
auto vec2<T>::lengthSquared() const -> T
{
    return x * x + y * y;
}

template<class T>
auto vec2<T>::dot(const vec2& _v) const -> T
{
    return x * _v.x + y * _v.y;
}

template<class T>
auto vec2<T>::cross(const vec2& _v) const -> T
{
    return x * _v.y - y * _v.x;
}

template<class T>
auto vec2<T>::max() const -> T
{
    return std::max(x, y);
}

template<class T>
auto vec2<T>::min() const -> T
{
    return std::min(x, y);
}

template<class T>
auto vec2<T>::max(const vec2& _v) const -> vec2
{
    return vec2(std::max(x, _v.x), std::max(y, _v.y));
}

template<class T>
auto vec2<T>::min(const vec2& _v) const -> vec2
{
    return vec2(std::min(x, _v.x), std::min(y, _v.y));
}

template struct vec2<float>;
template struct vec2<double>;

