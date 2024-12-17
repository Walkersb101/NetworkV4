#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include <range/v3/view/drop_last.hpp>
#include <range/v3/view/join.hpp>

#include "DataOut.hpp"

namespace IO
{
namespace timeSeries
{
class CSVOut : public dataOut
{
public:
  CSVOut(const std::filesystem::path& _timeDataPath,
         const std::filesystem::path& _bondDataPath);
  CSVOut(const std::filesystem::path& _timeDataPath,
         const std::filesystem::path& _bondDataPath,
         const std::string& _timeDataName,
         const std::string& _bondName);
  ~CSVOut();

public:
  void initFiles(const std::vector<std::string>& _dataHeader,
                 const std::vector<std::string>& _bondHeader) override;

  void writeTimeData(const std::vector<writeableTypes>& _data) override;
  void writeBondData(const std::vector<writeableTypes>& _data) override;

private:
  void writeHeader(const std::filesystem::path& _path,
                   const std::vector<std::string>& _header);
  void writeData(const std::filesystem::path& _path,
                 const std::vector<writeableTypes>& _data);

  std::string stripExtension(const std::string& _name);
};

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
  m_timeDataName = stripExtension(m_timeDataName) + config::IO::CSV::extension;
  writeHeader(m_timeDataPath / m_timeDataName, _dataHeader);

  m_bondName = stripExtension(m_bondName) + config::IO::CSV::extension;
  writeHeader(m_bondDataPath / m_bondName, _bondHeader);
}

void CSVOut::writeTimeData(const std::vector<writeableTypes>& _data)
{
  writeData(m_timeDataPath / m_timeDataName, _data);
}

void CSVOut::writeBondData(const std::vector<writeableTypes>& _data)
{
  writeData(m_bondDataPath / m_bondName, _data);
}

void CSVOut::writeHeader(const std::filesystem::path& _path,
                         const std::vector<std::string>& _header)
{
  std::ofstream file(_path);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _path.string());
  }
  for (const auto& x : _header | ranges::views::drop_last(1)) {
    file << x << ",";
  }
  file << _header.back() << "\n";
  file.close();
}

void CSVOut::writeData(const std::filesystem::path& _path,
                       const std::vector<writeableTypes>& _data)
{
  std::ofstream file(_path, std::ios::app);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _path.string());
  }
  file << std::setprecision(config::IO::CSV::precision);
  for (const auto& x : _data | ranges::views::drop_last(1)) {
    std::visit([&file](const auto& x) { file << x << ","; }, x);
  }
  std::visit([&file](const auto& x) { file << x << "\n"; }, _data.back());
  file.close();
}

std::string CSVOut::stripExtension(const std::string& _name)
{
  if (std::filesystem::path(_name).has_extension()) {
    return _name.substr(0, _name.find_last_of('.'));
  }
  return _name;
}

}  // namespace timeSeries
}  // namespace IO