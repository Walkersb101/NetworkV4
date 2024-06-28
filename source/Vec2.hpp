#pragma once

#include <cmath>

template <class T>
struct vec2 {
    T x, y;
    vec2();
    vec2(T _x, T _y);
    vec2(const vec2& _v);

public:
    auto operator=(const vec2& _v) -> vec2&;
    
    auto operator+(const T& _v) const -> const vec2;
    auto operator-(const T& _v) const -> const vec2;
    auto operator*(const T& _s) const -> const vec2;
    auto operator/(const T& _s) const -> const vec2;

    auto operator+(const vec2& _v) const -> const vec2;
    auto operator-(const vec2& _v) const -> const vec2;
    auto operator*(const vec2& _v) const -> const vec2;
    auto operator/(const vec2& _v) const -> const vec2;

    auto operator+=(const T& _v) -> vec2&;
    auto operator-=(const T& _v) -> vec2&;
    auto operator*=(const T& _s) -> vec2&;
    auto operator/=(const T& _s) -> vec2&;

    auto operator+=(const vec2& _v) -> vec2&;
    auto operator-=(const vec2& _v) -> vec2&;
    auto operator*=(const vec2& _v) -> vec2&;
    auto operator/=(const vec2& _v) -> vec2&;

    void set(const T& _x, const T& _y);
    void set(const vec2& _v);

    void normalize();
    auto length() const -> T;
    auto lengthSquared() const -> T;

    auto dot(const vec2& _v) const -> T;
    auto cross(const vec2& _v) const -> T;

    auto max() const -> T;
    auto min() const -> T;
    auto max(const vec2& _v) const -> vec2;
    auto min(const vec2& _v) const -> vec2;
};

using vec2f = vec2<float>;
using vec2d = vec2<double>;

#pragma omp declare reduction (+: vec2f: omp_out += omp_in) initializer(omp_priv = vec2f(0, 0))
#pragma omp declare reduction (+: vec2d: omp_out += omp_in) initializer(omp_priv = vec2d(0, 0))