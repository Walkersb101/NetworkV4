#include <fstream>
#include <iostream>
#include <stdexcept>

#include "NetworkOut.hpp"

IO::networkDumps::networkOut::networkOut() {}

IO::networkDumps::networkOut::networkOut(const std::filesystem::path& _path)
    : m_path(_path)
{
  if (!std::filesystem::exists(m_path.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_path.parent_path().string());
  }
}

IO::networkDumps::networkOut::~networkOut() {}

void IO::networkDumps::networkOut::save(const networkV4::network& _net,
                      const std::size_t _step,
                      const double _time,
                      const std::string& _type) const
{
}


