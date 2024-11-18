#pragma once

#include <cstdint>

#include "Core/BreakTypes.hpp"
#include "Core/Network.hpp"
#include "IO/NetworkOut.hpp"
#include "IO/DataOut.hpp"
#include "Misc/Enums.hpp"

namespace networkV4
{
class protocol
{
public:
  protocol();
  protocol(StrainType _strainType, std::unique_ptr<BreakTypes> _breakType);
  virtual ~protocol();

public:
  double getStrain(const network& _network) const;
  void strain(network& _network, double _step);

public:
  virtual void run(network& _network);

  virtual void initIO(const network& _network,
                      std::unique_ptr<dataOut>& _dataOut,
                      std::unique_ptr<networkOut>& _networkOut);

protected:
  StrainType m_strainType;
  BreakTypes* m_breakProtocol;
  dataOut* m_dataOut;
  networkOut* m_networkOut;
};

std::vector<double> forceMags(const network& _network);

}  // namespace networkV4
