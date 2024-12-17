#include "IO/NetworkOut.hpp"

namespace IO
{
namespace networkDumps
{
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