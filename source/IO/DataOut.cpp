#include "DataOut.hpp"

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

IO::timeSeries::dataOut::dataOut(const std::filesystem::path& _timeDataPath,
                 const std::filesystem::path& _bondDataPath,
                 const std::string& _timeDataName,
                 const std::string& _bondName)
    : m_timeDataPath(_timeDataPath)
    , m_bondDataPath(_bondDataPath)
    , m_timeDataName(_timeDataName)
    , m_bondName(_bondName)
{
  if (!std::filesystem::exists(m_timeDataPath.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_timeDataPath.parent_path().string());
  }
  if (!std::filesystem::exists(m_bondDataPath.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_bondDataPath.parent_path().string());
  }
}

IO::timeSeries::dataOut::~dataOut() {}

void IO::timeSeries::dataOut::initFiles(const std::vector<std::string>& _dataHeader,
                        const std::vector<std::string>& _bondHeader)
{
}


void IO::timeSeries::dataOut::dataOut::writeTimeData(const std::vector<writeableTypes>& _data) {}

void IO::timeSeries::dataOut::dataOut::writeBondData(const std::vector<writeableTypes>& _data) {}