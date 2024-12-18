#pragma once

#include <filesystem>
#include <memory>
#include <variant>

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Config.hpp"
#include "Misc/Enums.hpp"

namespace IO
{

using writeableTypes = std::variant<std::string, double, std::size_t>;

namespace timeSeries
{

class baseOut
{
public:
  baseOut();
  virtual ~baseOut() = 0;

public:
  virtual void write(const std::vector<writeableTypes>& _data) = 0;
};

class fileOut : public baseOut
{
public:
  fileOut(const std::filesystem::path& _dataPath,
          const std::string& _fileName,
          const std::vector<writeableTypes>& _header);

  virtual ~fileOut();

public:
  virtual void write(const std::vector<writeableTypes>& _data) = 0;

protected:
  std::filesystem::path m_filePath = "";
  size_t m_nCols = 0;
};

class NoOut : public baseOut
{
public:
  NoOut();
  ~NoOut();

public:
  void write(const std::vector<writeableTypes>& _data) override;
};

}  // namespace timeSeries
}  // namespace IO
