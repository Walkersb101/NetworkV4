#pragma once

#include <filesystem>
#include <variant>
#include <vector>

namespace IO
{

auto validFolder(const std::filesystem::path& _dirPath) -> bool;

class baseIO
{
protected:
  baseIO();
  virtual ~baseIO() = 0;
};

class folderIO : public baseIO
{
protected:
  folderIO(const std::filesystem::path& _dirPath);
  virtual ~folderIO() = 0;

protected:
    std::filesystem::path m_dirPath;
};

class fileIO : public baseIO
{
protected:
  fileIO(const std::filesystem::path& _filePath);
  virtual ~fileIO() = 0;

protected:
    std::filesystem::path m_filePath;
};
}  // namespace IO