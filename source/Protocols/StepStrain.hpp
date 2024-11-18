
#pragma once

#include <cstdint>

#include "Core/BreakTypes.hpp"
#include "Core/Network.hpp"
#include "IO/NetworkOut.hpp"
#include "IO/DataOut.hpp"
#include "Misc/Enums.hpp"

#include "Protocols/Protocol.hpp"

namespace networkV4
{

class stepStrain : public protocol
{
public:
  stepStrain();
  stepStrain(double _maxStrain,
             StrainType _strainType,
             std::unique_ptr<BreakTypes> _breakType,
             double _maxTime,
             double _esp,
             double _timeScale,
             double _stressScale);
  ~stepStrain();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  double globalStressAlongAxis(const network& _network) const;
  void checkLog(const network& _network);
  auto checkExit(const network& _network) -> bool;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   size_t _bondIndex,
                   BreakTypes* _breakProtocol) -> std::vector<writeableTypes>;

private:
  double m_maxStrain;

  double m_esp;
  double m_t;
  double m_maxTime;

  double m_timeScale;
  double m_timeLogPar;

  double m_stressScale;
  double m_stressLogPar;

  size_t m_logCount;
};

}  // namespace networkV4