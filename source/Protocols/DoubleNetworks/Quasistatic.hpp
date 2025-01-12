#pragma once

#include <cstdint>

#include "Integration/Integrators/Adaptive.hpp"
#include "Integration/Minimizers/AdaptiveHeunDecent.hpp"
#include "Integration/Minimizers/Fire2.hpp"
#include "Misc/Config.hpp"
#include "Misc/Roots.hpp"
#include "Protocols/Protocol.hpp"
#include "Protocols/deform.hpp"
#include "Protocols/protocolReader.hpp"

namespace networkV4
{
namespace protocols
{

class quasiStaticStrainDouble : public protocolBase
{
public:
  quasiStaticStrainDouble() = delete;
  quasiStaticStrainDouble(
      std::shared_ptr<deform::deformBase>& _deform,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
      std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
      const network& _network,
      double _maxStrain,
      double _rootTol = config::rootMethods::targetTol,
      integration::AdaptiveParams _params = integration::AdaptiveParams(),
      const minimisation::minimiserParams& _minParams =
          minimisation::minimiserParams(),
      bool _errorOnNotSingleBreak = false,
      double _maxStep = config::protocols::maxStep,
      bool _saveBreaks = false);
  ~quasiStaticStrainDouble();

public:
  void run(network& _network) override;

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
  auto evalStrain(network& _network,
                  double _targetStrain) -> std::tuple<double, size_t>;
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
                   std::size_t _breakCount)
      -> std::vector<IO::timeSeries::writeableTypes>;
  auto genBondData(const network& _network, const breakInfo& _bond)
      -> std::vector<IO::timeSeries::writeableTypes>;
  auto getCounts(const network& _network)
      -> std::tuple<std::size_t, std::size_t, std::size_t>;

private:
  auto breakData(const network& _network) -> std::tuple<double, size_t>;

private:
  auto xdoty(const std::vector<Utils::Math::vec2d>& _x,
             const std::vector<Utils::Math::vec2d>& _y) -> double;

private:
  double m_maxStrain;
  bool m_oneBreak = false;

  integration::AdaptiveParams m_params;
  minimisation::minimiserParams m_minParams;
  double m_rootTol;

  bool m_errorOnNotSingleBreak;

  double m_maxStep;

  double m_t = 0.0;
  std::size_t m_strainCount = 0;

  bool m_saveBreaks = false;
};

class quasiStaticStrainDoubleReader : public protocolReader
{
public:
  quasiStaticStrainDoubleReader() = default;
  ~quasiStaticStrainDoubleReader() = default;

public:
  auto read(const toml::value& _config,
            const network& _network,
            std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
            std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
            std::shared_ptr<IO::networkDumps::networkDump>& _networkOut)
      -> std::shared_ptr<protocolBase> override;
};

}  // namespace protocols

}  // namespace networkV4