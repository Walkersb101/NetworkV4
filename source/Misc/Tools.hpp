#pragma once

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <numeric>
#include <unordered_map>
#include <vector>

#include "Core/Bonds.hpp"
#include "Core/Network.hpp"
#include "EnumString.hpp"
#include "Enums.hpp"
#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace tools
{

template<typename T>
inline auto norm(const std::vector<T>& _vec) -> double
{
  double sum =
      std::accumulate(_vec.begin(),
                      _vec.end(),
                      0.0,
                      [](double _sum, const T& _v) { return _sum + _v * _v; });
  return std::sqrt(sum);
}

template<typename U>
inline auto norm(const std::vector<Utils::vec2<U>>& _vec) -> double
{
  const double sum = std::accumulate(_vec.begin(),
                               _vec.end(),
                               0.0,
                               [](double _sum, const Utils::vec2<U>& _v)
                               { return _sum + _v.normSquared(); });
  return std::sqrt(sum);
}

template<typename T>
inline auto maxLength(const std::vector<Utils::vec2<T>>& _vec) -> T
{
  // double max = 0.0;
  // for (size_t i = 0; i < _vec.size(); ++i) {
  //     if (_vec[i].length() > max)
  //         max = _vec[i].length();
  // }
  // return max;

  return std::accumulate(_vec.begin(),
                         _vec.end(),
                         0.0,
                         [](double _max, const Utils::vec2<T>& _v)
                         { return std::max(_max, _v.norm()); });
}


template<typename T>
inline auto maxAbsComponent(const std::vector<Utils::vec2<T>>& _vec) -> T
{
  return std::accumulate(_vec.begin(),
                         _vec.end(),
                         0.0,
                         [](double _max, const Utils::vec2<T>& _v)
                         { return std::max(_max, _v.abs().max()); });
}

inline auto xdoty(const std::vector<Utils::vec2d>& _vec1, const std::vector<Utils::vec2d>& _vec2) -> double
{
  return std::inner_product(_vec1.begin(),
                            _vec1.end(),
                            _vec2.begin(),
                            0.0,
                            std::plus<double>(),
                            [](const Utils::vec2d& _v1, const Utils::vec2d& _v2)
                            { return _v1.dot(_v2); });
}

template<typename T>
inline auto sign(const T& _val) -> T
{
  return _val < 0 ? -1 : 1;
}

template<typename T>
inline auto contains(const std::vector<T>& _vec, const T& _val) -> bool
{
  return std::find(_vec.begin(), _vec.end(), _val) != _vec.end();
}

inline void checkCanOpen(const std::filesystem::path& _filename)
{
  if (!std::filesystem::exists(_filename)) {
    throw std::runtime_error("File does not exist: " + _filename.string());
  }
  std::ifstream file(_filename, std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _filename.string());
  }
}

}  // namespace tools

namespace tensorData
{
inline void addBondCountByTypeHeader(
    std::vector<std::string>& _header,
    const std::vector<networkV4::bondType>& _types)
{
  _header.reserve(_header.size() + _types.size() + 1);
  _header.emplace_back("BondCount");
  for (const auto& type : _types) {
    const std::string bondType = enumString::bondType2Str.at(type);
    _header.emplace_back(bondType + "Count");
  }
}

inline void addStressByTypeHeader(
    std::vector<std::string>& _header,
    const std::vector<networkV4::bondType>& _types)
{
  const std::vector<std::string> header = {
      "Stressxx", "Stressxy", "Stressyx", "Stressyy"};
  _header.reserve(_header.size() + (_types.size() + 1) * 4);
  _header.insert(_header.end(), header.begin(), header.end());
  for (const auto& type : _types) {
    const std::string bondType = enumString::bondType2Str.at(type);
    for (const auto& col : header) {
      _header.emplace_back(bondType + col);
    }
  }
}

}  // namespace tensorData

}  // namespace networkV4