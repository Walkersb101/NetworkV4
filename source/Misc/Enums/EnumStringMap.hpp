#pragma once

#include <array>
#include <initializer_list>
#include <stdexcept>
#include <string>

#include <range/v3/algorithm/copy.hpp>

namespace Utils
{

template<typename EnumType, std::size_t N>
class EnumStringLookup
{
public:
  // Constructor accepting an initializer list
  EnumStringLookup(
      std::initializer_list<std::pair<EnumType, std::string>> enumStringList)
  {
    if (enumStringList.size() != N) {
      throw std::invalid_argument(
          "Initializer list size does not match the expected size.");
    }

    ranges::copy(enumStringList, enumString.begin());
  };

  // Lookup string by enum value
  inline std::string lookup(EnumType enumValue) const
  {
    for (const auto& pair : enumString) {
      if (pair.first == enumValue) {
        return pair.second;
      }
    }
    throw std::invalid_argument("Invalid enum value");
  };

  // Reverse lookup: Find enum by string value
  inline EnumType reverseLookup(const std::string& str) const
  {
    for (const auto& pair : enumString) {
      if (pair.second == str) {
        return pair.first;
      }
    }
    throw std::invalid_argument("Invalid string value");
  };

private:
  std::array<std::pair<EnumType, std::string>, N>
      enumString;  // Now an internal array
};

}; // namespace Utils
