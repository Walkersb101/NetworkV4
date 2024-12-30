#include <fstream>
#include <iostream>
#include <stdexcept>

#include "NetworkOut.hpp"

IO::networkDumps::networkDump::networkDump() {}

IO::networkDumps::networkDump::~networkDump() {}

void IO::networkDumps::networkDump::save(const networkV4::network& _net,
                                         const std::size_t _step,
                                         const double _time,
                                         const std::string& _type) const
{
}

IO::networkDumps::noNetworkOut::noNetworkOut() {}

IO::networkDumps::noNetworkOut::~noNetworkOut() {}

void IO::networkDumps::noNetworkOut::save(const networkV4::network& _net,
                                          const std::size_t _step,
                                          const double _time,
                                          const std::string& _type) const
{
}
