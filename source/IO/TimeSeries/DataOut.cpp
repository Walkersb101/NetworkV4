
#include "DataOut.hpp"

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

IO::timeSeries::baseOut::baseOut() {}

IO::timeSeries::baseOut::~baseOut() {}

IO::timeSeries::fileOut::fileOut(const std::filesystem::path& _dataPath,
                                 const std::string& _fileName,
                                 const std::vector<writeableTypes>& _dataHeader)
    : m_filePath(_dataPath / _fileName)
    , m_nCols(_dataHeader.size())
{
  if (!std::filesystem::exists(m_filePath.parent_path())) {
    std::string err =
        "Directory " + m_filePath.parent_path().string() + " does not exist.";
    throw std::runtime_error(err);
  }
}

IO::timeSeries::fileOut::~fileOut() {}

void IO::timeSeries::fileOut::write(const std::vector<writeableTypes>& _data) {}

IO::timeSeries::NoOut::NoOut() {}

IO::timeSeries::NoOut::~NoOut() {}

void IO::timeSeries::NoOut::write(const std::vector<writeableTypes>& _data) {}