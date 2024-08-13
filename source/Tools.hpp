#pragma once

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <bxzstr.hpp>

#include "Bonds.hpp"
#include "EnumString.hpp"
#include "Enums.hpp"
#include "Vec2.hpp"

namespace networkV4
{
namespace tools
{
template<typename T, typename U>
inline auto averageMaps(const std::unordered_map<T, U>& _map1,
                        const std::unordered_map<T, U>& _map2)
    -> std::unordered_map<T, U>
{
  std::unordered_map<T, U> averageMap;
  for (const auto& [key, val] : _map1) {
    if (_map2.find(key) == _map2.end()) {
      throw std::runtime_error("Key not found in _stresses2");
    } else {
      averageMap[key] = (val + _map2.at(key)) * 0.5;
    }
  }
  return averageMap;
}

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
inline auto norm(const std::vector<vec2<U>>& _vec) -> double
{
  double sum = std::accumulate(_vec.begin(),
                               _vec.end(),
                               0.0,
                               [](double _sum, const vec2<U>& _v)
                               { return _sum + _v.lengthSquared(); });
  return std::sqrt(sum);
}

template<typename T>
inline auto maxLength(const std::vector<vec2<T>>& _vec) -> T
{
  return std::accumulate(_vec.begin(),
                         _vec.end(),
                         0.0,
                         [](double _max, const vec2<T>& _v)
                         { return std::max(_max, _v.length()); });
}

inline auto uniqueBondTypes(const networkV4::bonds& _bonds)
    -> std::vector<networkV4::bondType>
{
  std::vector<networkV4::bondType> types;
  for (const auto& bond : _bonds) {
    if (std::find(types.begin(), types.end(), bond.type()) == types.end()) {
      types.push_back(bond.type());
    }
  }
  return types;
}

template<typename T>
inline auto maxAbsComponent(const std::vector<vec2<T>>& _vec) -> T
{
  return std::accumulate(_vec.begin(),
                         _vec.end(),
                         0.0,
                         [](double _max, const vec2<T>& _v)
                         { return std::max(_max, _v.abs().max()); });
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

} // namespace tensorData

#pragma omp declare reduction(vec_merge : std::vector<int> : omp_out.insert( \
        omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction( \
        vec_merge : std::vector<std::size_t> : omp_out.insert( \
                omp_out.end(), omp_in.begin(), omp_in.end()))
}  // namespace networkV4