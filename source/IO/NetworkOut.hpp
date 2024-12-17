#pragma once

#include <filesystem>

#include "Core/Network.hpp"

namespace IO
{
namespace networkDumps
{
class networkOut
{
public:
  networkOut();
  networkOut(const std::filesystem::path& _path);
  virtual ~networkOut();

public:
  virtual void save(const networkV4::network& _net,
                    const std::size_t _step,
                    const double _time,
                    const std::string& _type) const;

protected:
  std::filesystem::path m_path;
};

}  // namespace networkDumps
}  // namespace IO