#pragma once

#include <cstdint>
#include <stdexcept>
#include <string>
#include <optional>

#include "IO/Compression/zstd.hpp"
#include "Misc/Enums/EnumStringMap.hpp"

namespace IO
{
namespace Compression
{

enum class compressionTypes : uint8_t
{
  None,
  Zstd,
};

const Utils::EnumStringLookup<compressionTypes, 2> compressionExts = {
    { { compressionTypes::None, "" },
      { compressionTypes::Zstd, ".zst" } } };

class compression
{
public:
  compression(compressionTypes _compression, const std::optional<size_t> _level = std::nullopt)
      : m_type(_compression)
  {
    if (_level.has_value()) {
      m_level = _level.value();
      checkCompressionLevel();
    } else {
      setDefaultCompressionLevel();
    }
  }
  ~compression() {};

public:
  std::string compress(const std::string& _data) const
  {
    switch (m_type) {
      case compressionTypes::None:
        return _data;
      case compressionTypes::Zstd:
        return zstd::compress(_data, m_level);
      default:
        throw std::runtime_error("Unknown compression type");
    }
  }

private:
  void setDefaultCompressionLevel()
  {
    switch (m_type) {
      case compressionTypes::None:
        break;
      case compressionTypes::Zstd:
        m_level = 3;
        break;
      default:
        throw std::runtime_error("Unknown compression type");
    };
  }

  void checkCompressionLevel()
  {
    switch (m_type) {
      case compressionTypes::None:
        break;
      case compressionTypes::Zstd:
        if (m_level > 22) {
          throw std::runtime_error("Zstd compression level must be between 0 and 22"); //todo: should this throw an exception?
        }
        break;
      default:
        throw std::runtime_error("Unknown compression type");
    };
  }

protected:
  compressionTypes m_type;
  size_t m_level;
};

class decompression
{
public:
  decompression(compressionTypes _compression)
      : m_type(_compression)
  {
  }
  ~decompression() {};

public:
    std::string decompress(const std::string& _data) const
    {
        switch (m_type) {
        case compressionTypes::None:
            return _data;
        case compressionTypes::Zstd:
            return zstd::decompress(_data);
        default:
            throw std::runtime_error("Unknown compression type");
        }
    }

private:
    compressionTypes m_type;
};

inline auto detectCompression(const std::string& _data) -> compressionTypes
{
    if (zstd::isZstdCompressed(_data)) {
        return compressionTypes::Zstd;
    }
    return compressionTypes::None;
}

}  // namespace Compression
}  // namespace IO