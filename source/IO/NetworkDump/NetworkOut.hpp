#pragma once

#include <filesystem>

#include "Core/Network.hpp"

namespace IO
{
namespace networkDumps
{
class networkDump
{
public:
  networkDump();
  virtual ~networkDump() = 0;

public:
  virtual void save(const networkV4::network& _net,
                    const std::size_t _step,
                    const double _time,
                    const std::string& _type) const = 0;
};

class noNetworkOut : public networkDump
{
public:
  noNetworkOut();
  ~noNetworkOut() override;

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override;
};

}  // namespace networkDumps
}  // namespace IO