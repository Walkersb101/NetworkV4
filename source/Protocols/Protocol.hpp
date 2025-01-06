#pragma once

#include <cstdint>

#include "Core/Network.hpp"
#include "IO/NetworkDump/NetworkOut.hpp"
#include "IO/TimeSeries/DataOut.hpp"
#include "Misc/Enums.hpp"
#include "deform.hpp"

namespace networkV4
{
namespace protocols
{
class protocolBase
{
public:
  protocolBase() = delete;
  template<typename def>
  protocolBase(def _deform,
               std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
               std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
               const network& _network)
      : m_deform(std::make_unique<deform::deformBase>(_deform))
      , m_dataOut(_dataOut)
      , m_networkOut(_networkOut)
  {
  }
  virtual ~protocolBase() = 0;

public:
  virtual void run(network& _network) = 0;

protected:
  std::unique_ptr<deform::deformBase> m_deform;

  std::shared_ptr<IO::timeSeries::timeSeriesOut> m_dataOut;
  std::shared_ptr<IO::networkDumps::networkDump> m_networkOut;
};

std::vector<double> forceMags(const network& _network);

}  // namespace protocols
}  // namespace networkV4
