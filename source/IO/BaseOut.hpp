#pragma once

#include <filesystem>
#include <variant>
#include <vector>

namespace IO
{

auto validFolder(const std::filesystem::path& _dirPath) -> bool;

class baseOut
{
protected:
  baseOut();
  virtual ~baseOut() = 0;
};

class folderOut : public baseOut
{
protected:
  folderOut(const std::filesystem::path& _dirPath);
  virtual ~folderOut() = 0;

protected:
    std::filesystem::path m_dirPath;
};

class fileOut : public baseOut
{
protected:
  fileOut(const std::filesystem::path& _filePath);
  virtual ~fileOut() = 0;

protected:
    std::filesystem::path m_filePath;
};
}  // namespace IO