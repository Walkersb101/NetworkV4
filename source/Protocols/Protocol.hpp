#pragma once

#include <cstdint>

#include "Core/Network.hpp"
#include "IO/DataOut.hpp"
#include "IO/NetworkOut.hpp"
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
                      std::unique_ptr<IO::timeSeries::dataOut>& _dataOut,
                      std::unique_ptr<IO::networkDumps::networkOut>& _networkOut);

protected:
  StrainType m_strainType;
  IO::timeSeries::dataOut* m_dataOut;
  IO::networkDumps::networkOut* m_networkOut;
};

std::vector<double> forceMags(const network& _network);

}  // namespace networkV4
