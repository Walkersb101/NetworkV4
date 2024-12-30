
#include "DataOut.hpp"

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

IO::timeSeries::timeSeriesOut::timeSeriesOut() {}

IO::timeSeries::timeSeriesOut::~timeSeriesOut() {}

void IO::timeSeries::timeSeriesOut::write(const std::vector<writeableTypes>& _data) {}

IO::timeSeries::NoOut::NoOut() {}

IO::timeSeries::NoOut::~NoOut() {}

void IO::timeSeries::NoOut::write(const std::vector<writeableTypes>& _data) {}