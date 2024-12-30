#pragma once

#include <cstdint>

#include "Core/Network.hpp"
#include "IO/NetworkDump/NetworkOut.hpp"
#include "IO/TimeSeries/DataOut.hpp"
#include "Misc/Enums.hpp"

namespace networkV4
{
class protocol
{
public:
  protocol();
  protocol(StrainType _strainType);
  virtual ~protocol();

public:
  double getStrain(const network& _network) const;
  void strain(network& _network, double _step);

public:
  virtual void run(network& _network);

  virtual void initIO(const network& _network,
                      std::unique_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
                      std::unique_ptr<IO::networkDumps::networkDump>& _networkOut);

protected:
  StrainType m_strainType;
  IO::timeSeries::timeSeriesOut* m_dataOut;
  IO::networkDumps::networkDump* m_networkOut;
};

std::vector<double> forceMags(const network& _network);

}  // namespace networkV4
