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

class quasiStaticStrain : public protocol
{
public:
  quasiStaticStrain();
  quasiStaticStrain(double _maxStrain,
                    StrainType _strainType,
                    std::unique_ptr<BreakTypes>& _breakType);
  quasiStaticStrain(double _maxStrain,
                    StrainType _strainType,
                    std::unique_ptr<BreakTypes>& _breakType,
                    double _esp,
                    double _tol,
                    bool _errorOnNotSingleBreak,
                    double _maxStep,
                    bool _single = false);
  ~quasiStaticStrain();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  void evalStrain(network& _network,
                  double _step,
                  double& _maxVal,
                  std::size_t& _count,
                  bool _save = false);
  auto converge(network& _baseNetwork,
                double& _a,
                double& _b,
                double& _fa,
                double& _fb,
                size_t& _breakCountB,
                double _tol,
                bool _earlyExit = false) -> bool;
  auto findSingleBreak(network& _network) -> size_t;

  auto relaxBreak(network& _network) -> std::size_t;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   size_t _bondIndex) -> std::vector<writeableTypes>;

private:
  double m_maxStrain;
  bool m_single;

  double m_esp;
  double m_tol;
  bool m_errorOnNotSingleBreak;
  double m_maxStep;

  double m_t;

  std::size_t m_strainCount;

  bool m_saveBreaks;
};

}  // namespace networkV4
