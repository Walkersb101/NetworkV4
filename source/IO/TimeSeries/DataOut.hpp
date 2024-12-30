#pragma once

#include <filesystem>
#include <memory>
#include <variant>

#include "IO/BaseIO.hpp"

namespace IO
{
namespace timeSeries
{

using writeableTypes = std::variant<std::string, double, std::size_t>;

class timeSeriesOut
{
public:
  timeSeriesOut();
  virtual ~timeSeriesOut() = 0;

public:
  virtual void write(const std::vector<writeableTypes>& _data) = 0;
};

class NoOut : public timeSeriesOut
{
public:
  NoOut();
  ~NoOut();

public:
  void write(const std::vector<writeableTypes>& _data) override;
};

}  // namespace timeSeries
}  // namespace IO
