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
  virtual ~networkOut() = 0;

public:
  virtual void save(const networkV4::network& _net,
                    const std::size_t _step,
                    const double _time,
                    const std::string& _type) const = 0;

protected:
  std::filesystem::path m_path;
};

class noNetworkOut : public networkOut
{
public:
  noNetworkOut() {};
  ~noNetworkOut() override {};

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override {};
};

}  // namespace networkDumps
}  // namespace IO