#pragma once

#include <optional>
#include <span>
#include <variant>
#include <vector>

#include "Core/Network.hpp"
#include "Misc/Math/Vector.hpp"

namespace Utils
{

template<typename T, std::size_t N>
auto spanView(const std::vector<Math::Vector<T, N>>& _vec) -> std::span<const T>
{
  return std::span<const T>(reinterpret_cast<const T*>(_vec.data()),
                            _vec.size() * N);
}

template<typename T, std::size_t N>
auto spanView(std::vector<Math::Vector<T, N>>& _vec) -> std::span<T>
{
  return std::span<T>(reinterpret_cast<T*>(_vec.data()), _vec.size() * N);
}

template<typename T>
inline auto sign(const T& _val) -> T
{
  return _val < 0 ? -1 : 1;
}

template<typename T, std::size_t N>
inline auto rms(const std::vector<Math::Vector<T, N>>& _vec) -> T
{
  T sum = 0.0;
  for (const auto& val : _vec) {
    sum += val.norm2();
  }
  return std::sqrt(sum / _vec.size());
}

template<typename T, std::size_t N>
inline auto maxComp(const std::vector<Math::Vector<T, N>>& _vec) -> T
{
  T max = 0.0;
  for (const auto& val : _vec) {
    max = std::max(max, val.abs().max());
  }
  return max;
}

inline auto brokenBonds(const networkV4::network& _net) -> std::size_t
{
  auto brokenFlag = _net.getTags().get("broken");
  const auto& tags = _net.getBonds().getTags();
  return std::count_if(tags.begin(),
                       tags.end(),
                       [&brokenFlag](const auto& tag)
                       { return Utils::Tags::hasTag(tag, brokenFlag); });
}

namespace Output
{

enum class stepScale
{
  Linear,
  Log
};

template<typename T>
struct linearStepInfo
{
  linearStepInfo(T _initial, T _step)
      : startValue(_initial)
      , nextLog(_initial)
      , stepSize(_step)
  {
    assert(_step > 0);
  }
  T startValue;
  T nextLog;
  T stepSize;

  void update() { nextLog += stepSize; }
  void reset() { nextLog = startValue; }
};

template<typename T, typename U>
struct logStepInfo
{
  logStepInfo(T _initial, U _step)
      : startValue(_initial)
      , nextLog(_initial)
      , stepSize(_step)
  {
    assert(_step > 0);
    assert(_initial > 0);
    update();
  }
  T startValue;
  T nextLog;
  U stepSize;

  void update() { nextLog *= stepSize; }
  void reset() { nextLog = startValue; }
};

template<typename T, typename U>
using stepInfo = std::variant<linearStepInfo<T>, logStepInfo<T, U>>;

template<typename T, typename U>
constexpr bool operator>=(T _a, stepInfo<T, U>& _b)
{
  return std::visit(
      [_a](auto& info)
      {
        if (_a >= info.nextLog) {
          info.update();
          return true;
        }
        return false;
      },
      _b);
}

template<typename T, typename U>
constexpr bool operator>=(T _a, std::optional<stepInfo<T, U>>& _b)
{
  return _b.has_value() && _a >= _b.value();
}

}  // namespace Output

}  // namespace Utils