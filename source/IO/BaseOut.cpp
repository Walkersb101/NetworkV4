#include <fstream>

#include "IO/BaseOut.hpp"

auto IO::validFolder(const std::filesystem::path& _dirPath) -> bool
{
  if (!std::filesystem::is_directory(_dirPath)) {
    return false;
  }
  auto perms = std::filesystem::status(_dirPath).permissions();
  return (std::filesystem::perms::owner_read & perms)
      == std::filesystem::perms::owner_read
      && (std::filesystem::perms::owner_write & perms)
      == std::filesystem::perms::owner_write;
}

IO::baseOut::baseOut() {}

IO::baseOut::~baseOut() {}

IO::folderOut::folderOut(const std::filesystem::path& _dirPath)
    : baseOut()
    , m_dirPath(_dirPath)
{
  if (!validFolder(_dirPath)) {
    std::string err = "Directory " + _dirPath.string() + " is not valid.";
    throw std::runtime_error(err);
  }
}

IO::folderOut::~folderOut() {}

IO::fileOut::fileOut(const std::filesystem::path& _filePath)
    : baseOut()
    , m_filePath(_filePath)
{
  // check if the parent directory exists
  if (!std::filesystem::exists(_filePath.parent_path())) {
    std::string err =
        "Directory " + _filePath.parent_path().string() + " does not exist.";
    throw std::runtime_error(err);
  }
  // check if the file can be opened in write mode
  std::ofstream ofs(_filePath, std::ios::out);
  if (!ofs) {
    std::string err =
        "File " + _filePath.string() + " cannot be opened in write mode.";
    throw std::runtime_error(err);
  }
  ofs.close();
}
