#pragma once

#include <algorithm>
#include <string>
#include <unordered_map>
#include <vector>

#include "Misc/Tags/TagMap.hpp"

namespace Utils
{
namespace Tags
{

using tagFlags = std::bitset<NUM_TAGS>;
using tagStorage = std::vector<tagFlags>;

inline auto hasTag(const tagFlags& _tags, const tagFlags& _tag) -> bool
{
  return (_tags & _tag) == _tag;
}

inline auto hasTagAny(const tagFlags& _tags, const tagFlags& _tag) -> bool
{
  return (_tags & _tag).any();
}

}  // namespace Tags
}  // namespace Utils