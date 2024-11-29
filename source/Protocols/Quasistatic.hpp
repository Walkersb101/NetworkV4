#pragma once

#include <cstdint>

#include "Core/BreakTypes.hpp"
#include "Core/Network.hpp"
#include "IO/DataOut.hpp"
#include "IO/NetworkOut.hpp"
#include "Integration/Minimization.hpp"
#include "Misc/Enums.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"
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
                    bool _saveBreaks = false);
  ~quasiStaticStrain();

public:
  void run(network& _network) override;
  void initIO(const network& _network,
              std::unique_ptr<dataOut>& _dataOut,
              std::unique_ptr<networkOut>& _networkOut);

private:
  enum class SingleBreakReason : std::uint8_t
  {
    NoBreaksInStep,
    BreakAtLowerBound,
    DidNotConverge,
    MaxStrainReached,
    Complete
  };

private:
  void evalStrain(network& _network,
                  double _step,
                  double& _targetStrain,
                  std::size_t& _count);
  auto converge(network& _baseNetwork,
                double& _a,
                double& _b,
                double& _fa,
                double& _fb,
                size_t& _breakCountB,
                double _tol,
                bool _earlyExit = false) -> bool;
  auto findSingleBreak(network& _network)
      -> std::tuple<SingleBreakReason, size_t>;

  auto relaxBreak(network& _network) -> std::size_t;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount) -> std::vector<writeableTypes>;
  auto genBondData(const network& _network,
                   size_t _bondIndex) -> std::vector<writeableTypes>;

private:
  double m_maxStrain;
  bool m_oneBreak = false;

  double m_esp = config::integrators::adaptiveIntegrator::esp;
  double m_rootTol = config::rootMethods::targetTol;
  double m_forceTol = config::integrators::miminizer::tol;
  bool m_errorOnNotSingleBreak =
      config::protocols::quasiStaticStrain::errorOnNotSingleBreak;

  double m_maxStep = 0.0;
  double m_strainGuessScale =
      config::protocols::quasiStaticStrain::strainGuessScale;

  double m_t = 0.0;

  std::size_t m_strainCount = 0;

  bool m_saveBreaks = false;
};

}  // namespace networkV4