#pragma once

#include <string>

#include <zstd.h>

namespace IO
{
namespace Compression
{
namespace zstd
{
const std::string version = ZSTD_versionString();

inline std::string compress(const std::string& _data,
                            const int _compressionLevel = 3)
{
  const size_t compressedSize = ZSTD_compressBound(_data.size());
  std::string compressed(compressedSize, '\0');
  const size_t compressedSizeActual = ZSTD_compress(&compressed[0],
                                                    compressedSize,
                                                    _data.c_str(),
                                                    _data.size(),
                                                    _compressionLevel);
  compressed.resize(compressedSizeActual);
  compressed.shrink_to_fit();
  return compressed;
}

inline std::string decompress(const std::string& _data)
{
  const size_t decompressedSize =
      ZSTD_getFrameContentSize(_data.c_str(), _data.size());
  std::string decompressed(decompressedSize, '\0');
  const size_t decompressedSizeActual = ZSTD_decompress(
      &decompressed[0], decompressedSize, _data.c_str(), _data.size());
  decompressed.resize(decompressedSizeActual);
  decompressed.shrink_to_fit();
  return decompressed;
}

inline auto isZstdCompressed(const std::string& _data)
{
  const unsigned char b0 =
      *reinterpret_cast<const unsigned char*>(_data.data() + 0);
  const unsigned char b1 =
      *reinterpret_cast<const unsigned char*>(_data.data() + 1);
  const unsigned char b2 =
      *reinterpret_cast<const unsigned char*>(_data.data() + 2);
  const unsigned char b3 =
      *reinterpret_cast<const unsigned char*>(_data.data() + 3);
  return (b0 == 0x28 && b1 == 0xB5 && b2 == 0x2F && b3 == 0xFD);
}
}  // namespace zstd
}  // namespace Compression
}  // namespace IO