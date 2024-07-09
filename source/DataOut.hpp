#pragma once

#include <filesystem>
#include <memory>
#include <variant>

#include <toml.hpp>

#include "Network.hpp"
#include "Enums.hpp"

using writeableTypes = std::variant<std::string, double, std::size_t>;


class dataOut
{
public:
  dataOut(const std::filesystem::path& _timeDataPath,
          const std::filesystem::path& _bondDataPath,
          const std::string& _timeDataName,
          const std::string& _bondName);
  dataOut(const std::filesystem::path& _timeDataPath,
          const std::filesystem::path& _bondDataPath);
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
};
