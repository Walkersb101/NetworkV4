#include "DataOut.hpp"

#include <toml.hpp>

#include "Config.hpp"
#include "Network.hpp"
#include "Protocol.hpp"

dataOut::dataOut(const std::filesystem::path& _timeDataPath,
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

dataOut::dataOut(const std::filesystem::path& _timeDataPath,
                 const std::filesystem::path& _bondDataPath)
    : dataOut(_timeDataPath,
              _bondDataPath,
              config::IO::timeDataName,
              config::IO::bondDataName)
{
}

void dataOut::initFiles(const std::vector<std::string>& _dataHeader,
                        const std::vector<std::string>& _bondHeader)
{
}

dataOut::~dataOut() {}

void dataOut::writeTimeData(const std::vector<writeableTypes>& _data) {}

void dataOut::writeBondData(const std::vector<writeableTypes>& _data) {}

CSVOut::CSVOut(const std::filesystem::path& _timeDataPath,
               const std::filesystem::path& _bondDataPath)
    : dataOut(_timeDataPath, _bondDataPath)
{
}

CSVOut::CSVOut(const std::filesystem::path& _timeDataPath,
               const std::filesystem::path& _bondDataPath,
               const std::string& _timeDataName,
               const std::string& _bondName)
    : dataOut(_timeDataPath, _bondDataPath, _timeDataName, _bondName)
{
}

CSVOut::~CSVOut() {}

void CSVOut::initFiles(const std::vector<std::string>& _dataHeader,
                       const std::vector<std::string>& _bondHeader)
{
    std::string timeDataName = m_timeDataName;
    if (std::filesystem::path(m_timeDataName).has_extension()) {
        timeDataName = m_timeDataName.substr(0, m_timeDataName.find_last_of('.'));
    }
    m_timeDataName = timeDataName + ".csv";
    std::ofstream timeFile(m_timeDataPath / m_timeDataName);
    if (!timeFile.is_open()) {
        throw std::runtime_error("Cannot open file: " + m_timeDataPath.string());
    }
    for (std::size_t i = 0; i < _dataHeader.size() - 1; ++i) {
        timeFile << _dataHeader[i] << ",";
    }
    timeFile << _dataHeader.back() << "\n";
    timeFile.close();

    std::string bondDataName = m_bondName;
    if (std::filesystem::path(m_bondName).has_extension()) {
        bondDataName = m_bondName.substr(0, m_bondName.find_last_of('.'));
    }
    m_bondName = bondDataName + ".csv";
    std::ofstream bondFile(m_bondDataPath / m_bondName);
    if (!bondFile.is_open()) {
        throw std::runtime_error("Cannot open file: " + m_bondDataPath.string());
    }
    for (std::size_t i = 0; i < _bondHeader.size() - 1; ++i) {
        bondFile << _bondHeader[i] << ",";
    }
    bondFile << _bondHeader.back() << "\n";
    bondFile.close();
}

void CSVOut::writeTimeData(const std::vector<writeableTypes>& _data)
{
    std::ofstream timeFile(m_timeDataPath / m_timeDataName, std::ios::app);
    if (!timeFile.is_open()) {
        throw std::runtime_error("Cannot open file: " + m_timeDataPath.string());
    }
    timeFile << std::setprecision(config::IO::precision);
    for (std::size_t i = 0; i < _data.size() - 1; ++i) {
        std::visit([&timeFile](const auto &x) { timeFile << x << ","; }, _data[i]);
    }
    std::visit([&timeFile](const auto &x) { timeFile << x << "\n"; }, _data.back());
    timeFile.close();
}

void CSVOut::writeBondData(const std::vector<writeableTypes>& _data)
{
    std::ofstream bondFile(m_bondDataPath / m_bondName, std::ios::app);
    if (!bondFile.is_open()) {
        throw std::runtime_error("Cannot open file: " + m_bondDataPath.string());
    }
    bondFile << std::setprecision(config::IO::precision);
    for (std::size_t i = 0; i < _data.size() - 1; ++i) {
        std::visit([&bondFile](const auto &x) { bondFile << x << ","; }, _data[i]);
    }
    std::visit([&bondFile](const auto &x) { bondFile << x << "\n"; }, _data.back());
    bondFile.close();
}