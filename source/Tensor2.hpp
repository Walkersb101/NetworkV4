#pragma once

#include "Vec2.hpp"

template<class T>
struct tensor2{
public:
    T xx;
    T xy;
    T yx;
    T yy;
    tensor2();
    tensor2(T _xx, T _xy, T _yx, T _yy);
    tensor2(const tensor2& _t);

    auto operator=(const tensor2& _t) -> tensor2&;

    auto operator+(const T& _v) const -> const tensor2;
    auto operator-(const T& _v) const -> const tensor2;
    auto operator*(const T& _s) const -> const tensor2;
    auto operator/(const T& _s) const -> const tensor2;

    auto operator+(const tensor2& _v) const -> const tensor2;
    auto operator-(const tensor2& _v) const -> const tensor2;
    auto operator*(const tensor2& _v) const -> const tensor2;
    auto operator/(const tensor2& _v) const -> const tensor2;

    auto operator+=(const T& _v) -> tensor2&;
    auto operator-=(const T& _v) -> tensor2&;
    auto operator*=(const T& _s) -> tensor2&;
    auto operator/=(const T& _s) -> tensor2&;

    auto operator+=(const tensor2& _v) -> tensor2&;
    auto operator-=(const tensor2& _v) -> tensor2&;
    auto operator*=(const tensor2& _v) -> tensor2&;
    auto operator/=(const tensor2& _v) -> tensor2&;

    void set(const T& _xx, const T& _xy, const T& _yx, const T& _yy);
    void set(const tensor2& _t);
    void set(const T& _v);
};

template<typename T>
tensor2<T> outer(const vec2<T>& _v1, const vec2<T>& _v2);

using tensor2f = tensor2<float>;
using tensor2d = tensor2<double>;

#pragma omp declare reduction (+: vec2f: omp_out += omp_in) initializer(omp_priv = vec2f(0, 0))
#pragma omp declare reduction (+: vec2d: omp_out += omp_in) initializer(omp_priv = vec2d(0, 0))