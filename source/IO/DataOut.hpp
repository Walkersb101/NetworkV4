#pragma once

#include <filesystem>
#include <memory>
#include <variant>

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Enums.hpp"
#include "Misc/Config.hpp"

namespace IO
{

using writeableTypes = std::variant<std::string, double, std::size_t>;

namespace timeSeries
{ 
class dataOut
{
public:
  dataOut(
      const std::filesystem::path& _timeDataPath,
      const std::filesystem::path& _bondDataPath,
      const std::string& _timeDataName = config::IO::timeSeries::timeDataName,
      const std::string& _bondName = config::IO::timeSeries::bondDataName);

  virtual ~dataOut();

protected:
  std::filesystem::path m_timeDataPath;
  std::filesystem::path m_bondDataPath;
  std::string m_timeDataName;
  std::string m_bondName;

public:
  virtual void initFiles(const std::vector<std::string>& _dataHeader,
                         const std::vector<std::string>& _bondHeader);

  virtual void writeTimeData(const std::vector<writeableTypes>& _data);
  virtual void writeBondData(const std::vector<writeableTypes>& _data);
};

}  // namespace timeSeries
}  // namespace IO
