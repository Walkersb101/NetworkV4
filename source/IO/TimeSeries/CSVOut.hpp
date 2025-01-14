#pragma once

#include <filesystem>
#include <fstream>
#include <vector>

#include <range/v3/view/drop_last.hpp>
#include <range/v3/view/join.hpp>

#include "IO/BaseIO.hpp"
#include "IO/TimeSeries/DataOut.hpp"

namespace IO
{
namespace timeSeries
{
class CSVOut
    : public fileIO
    , public timeSeriesOut
{
public:
  CSVOut(const std::filesystem::path& _filePath)
      : fileIO(_filePath)
      , timeSeriesOut()
  {
    if (_filePath.extension() != ".csv") {
      throw std::invalid_argument("File extension must be .csv");
    }
    std::ofstream file(m_filePath, std::ios::out);
    if (!file.is_open()) {
      throw std::runtime_error(
          "Cannot open file: "
          + m_filePath.string());  // TODO: maybe this shouldnt just throw
    }
    file.close();
  }
  ~CSVOut() override = default;

public:
  void write(const std::vector<writeableTypes>& _data) override
  {
    writeRow(_data, std::ios::out | std::ios::app);
  }

private:
  void writeRow(const std::vector<writeableTypes>& _row,
                std::ios::openmode _mode)
  {
    std::ofstream file(m_filePath, _mode);
    if (!file.is_open()) {
      throw std::runtime_error(
          "Cannot open file: "
          + m_filePath.string());  // TODO: maybe this shouldnt just throw
    }
    file << std::setprecision(config::IO::CSV::precision);
    for (const auto& x : _row | ranges::views::drop_last(1)) {
      std::visit([&file](const auto& x) { file << x << ","; }, x);
    }
    std::visit([&file](const auto& x) { file << x << "\n"; }, _row.back());
    file.close();
  }
};

}  // namespace timeSeries
}  // namespace IO