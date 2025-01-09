#pragma once

#include "Misc/Math/Vector.hpp"

namespace Utils
{
template<class T>
struct vec2 : public Math::Vector<T, 2>
{
    using Base = Math::Vector<T, 2>;
    using Base::Base;
    vec2() = default;
    vec2(vec2 const&) = default;
    vec2& operator=(vec2 const&) = default;

    union {
        struct { T x, y; };
        typename Base::array_type data;
    };

    vec2(T x_, T y_) : data{x_, y_} {}
    vec2(const Base& base) : Base(base) {}

    auto dot(const vec2& _vec) const -> T
    {
        return (*this) * _vec;
    }
};

using vec2f = vec2<float>;
using vec2d = vec2<double>;
}  // namespace Utils