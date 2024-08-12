#pragma once

#include <algorithm>
#include <filesystem>
#include <numeric>
#include <unordered_map>
#include <vector>

#include <bxzstr.hpp>

#include "Bonds.hpp"
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

}  // namespace tools

namespace enum2str
{
inline auto bondTypeString(const networkV4::bondType& _type) -> std::string
{
  switch (_type) {
    case networkV4::bondType::single:
      return "Single";
    case networkV4::bondType::matrix:
      return "Matrix";
    case networkV4::bondType::sacrificial:
      return "Sacrificial";
    default:
      throw std::runtime_error("Unknown bond type");
  }
}

inline auto strainTypeString(const networkV4::StrainType& _type) -> std::string
{
  switch (_type) {
    case networkV4::StrainType::Shear:
      return "Shear";
    case networkV4::StrainType::Elongation:
      return "Elongation";
    default:
      throw std::runtime_error("Unknown strain type");
  }
}

inline void addBondCountByTypeHeader(
    std::vector<std::string>& _header,
    const std::vector<networkV4::bondType>& _types)
{
  _header.reserve(_header.size() + _types.size() + 1);
  _header.push_back("BondCount");
  for (const auto& type : _types) {
    _header.push_back(enum2str::bondTypeString(type) + "Count");
  }
}

inline void addStressByTypeHeader(
    std::vector<std::string>& _header,
    const std::vector<networkV4::bondType>& _types)
{
  _header.reserve(_header.size() + (_types.size() + 1) * 4);
  _header.insert(_header.end(),
                 {"Stressxx", "Stressxy", "Stressyx", "Stressyy"});
  for (const auto& type : _types) {
    _header.insert(_header.end(),
                   {enum2str::bondTypeString(type) + "Stressxx",
                    enum2str::bondTypeString(type) + "Stressxy",
                    enum2str::bondTypeString(type) + "Stressyx",
                    enum2str::bondTypeString(type) + "Stressyy"});
  }
}

inline auto compressionExt(const bxz::Compression _type) -> std::string
{
  // z, bz2, lzma, zstd, plaintext
  switch (_type) {
    case (bxz::Compression::z):
      return ".zz";
    case (bxz::Compression::bz2):
      return ".bz2";
    case (bxz::Compression::lzma):
      return ".lzma";
    case (bxz::Compression::zstd):
      return ".zst";
    case (bxz::Compression::plaintext):
      return "";
    default:
      throw std::runtime_error("Unknown compression type");
  }
}

inline auto outputExt(const networkOutType _type) -> std::string
{
  // z, bz2, lzma, zstd, plaintext
  switch (_type) {
    case (networkOutType::BinV2):
      return ".v2.bin";
    case (networkOutType::None):
      return "";
    default:
      throw std::runtime_error("Unknown compression type");
  }
}

}  // namespace enum2str

#pragma omp declare reduction (vec_merge : std::vector<int> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
#pragma omp declare reduction (vec_merge : std::vector<std::size_t> : omp_out.insert(omp_out.end(), omp_in.begin(), omp_in.end()))
}  // namespace networkV4