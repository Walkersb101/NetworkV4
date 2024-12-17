#pragma once

#include <cstdint>
#include <string>

#include "IO/Compression/zstd.hpp"
#include "Misc/Enums/EnumStringMap.hpp"

namespace IO
{
namespace Compression
{

enum class Compression : uint8_t
{
  None,
  Zstd,
};

const Utils::EnumStringLookup<Compression, 2> compressionExt = {
    { { Compression::None, "" }, { Compression::Zstd, ".zst" } } };

inline std::string compress(const std::string& _data, Compression _compression)

}  // namespace Compression
}  // namespace IO