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
class CSVOut : public fileOut
{
public:
  CSVOut(const std::filesystem::path& _dataPath,
         const std::string& _fileName,
         const std::vector<writeableTypes>& _header)
      : fileOut(_dataPath, stripExtension(_fileName) + ".csv", _header)
  {
    writeRow(_header, std::ios::out | std::ios::trunc);
  }
  ~CSVOut() override;

public:
  void write(const std::vector<writeableTypes>& _data) override
  {
    writeRow(_data, std::ios::out | std::ios::app);
  }

private:
  inline void writeRow(const std::vector<writeableTypes>& _row,
                       std::ios::openmode _mode)
  {
    std::ofstream file(m_filePath, _mode);
    if (!file.is_open()) {
      throw std::runtime_error(
          "Cannot open file: "
          + m_filePath.string());  // TODO: maybe this shouldnt just throw
    }
    for (const auto& x : _row | ranges::views::drop_last(1)) {
      std::visit([&file](const auto& x) { file << x << ","; }, x);
    }
    std::visit([&file](const auto& x) { file << x << "\n"; }, _row.back());
    file.close();
  }

  std::string stripExtension(const std::string& _name)
  {
    if (std::filesystem::path(_name).has_extension()) {
      return _name.substr(0, _name.find_last_of('.'));
    }
    return _name;
  }
};

}  // namespace timeSeries
}  // namespace IO