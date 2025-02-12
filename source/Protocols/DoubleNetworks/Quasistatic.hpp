#pragma once

#include <cstdint>

#include <tl/expected.hpp>

#include "Integration/Integrators/Adaptive.hpp"
#include "Integration/Minimizers/AdaptiveHeunDecent.hpp"
#include "Integration/Minimizers/Fire2.hpp"
#include "Integration/Minimizers/SD.hpp"
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
  struct networkSavePoints
  {
    networkSavePoints()
        : strainStep(std::nullopt)
        , breakCount(std::nullopt)
        , time(std::nullopt)
        , strain(std::nullopt)
    {
    }

    std::optional<Utils::Output::stepInfo<size_t, double>> strainStep =
        std::nullopt;
    std::optional<Utils::Output::stepInfo<size_t, double>> breakCount =
        std::nullopt;
    std::optional<Utils::Output::stepInfo<double, double>> time = std::nullopt;
    std::optional<Utils::Output::stepInfo<double, double>> strain =
        std::nullopt;
  };

  enum class nextBreakState : std::uint8_t
  {
    FoundSingleBreak,
    FoundMultipleBreaks,
    MaxStrainReached,
  };

  class relaxBreak : public minimisation::minimiserBase
  {
  public:
    relaxBreak() = delete;
    relaxBreak(quasiStaticStrainDouble& _protocol)
        : minimiserBase(_protocol.m_minParams)
        , m_params(_protocol.m_params)
        , m_protocol(_protocol)
    {
    }

  public:
    void minimise(network& _network) override;

  private:
    auto hybridStep(network& _network, auto& _stepper)
        -> tl::expected<double, lineSearch::lineSearchState>;

  private:
    integration::AdaptiveParams m_params;
    quasiStaticStrainDouble& m_protocol;
  };

private:
  double m_maxStrain;
  bool m_oneBreak = false;

  integration::AdaptiveParams m_params;
  minimisation::minimiserParams m_minParams;
  double m_rootTol;

  bool m_errorOnNotSingleBreak;

  double m_maxStep;

  std::size_t m_strainCount = 0;

  networkSavePoints m_savePoints;
  relaxBreak m_breakMinimiser;

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
      networkSavePoints _savePoints = networkSavePoints());
  ~quasiStaticStrainDouble();

public:
  void run(network& _network) override;

private:
  auto evalStrain(const network& _network, double _targetStrain) -> network;
  auto converge(network& _network,
                double _a,
                double _b,
                double _fa,
                double _fb,
                double _tol) -> roots::rootState;

  auto findNextBreak(network& _network) -> nextBreakState;

  void relaxBreak(network& _network, size_t startBreaks);

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount,
                   double _t = 0.0)
      -> std::vector<IO::timeSeries::writeableTypes>;
  auto genBondData(const network& _network, const breakInfo& _bond, double _t)
      -> std::vector<IO::timeSeries::writeableTypes>;
  void logData(network& _network,
               const std::string& _reason,
               double _breakcount,
               double _t,
               bool _writeDump);

  auto getCounts(const network& _network)
      -> std::tuple<std::size_t, std::size_t, std::size_t>;
  auto checkIfNeedToSave(const network& _network, double _t = 0.0)
      -> std::optional<std::string>;
  auto processBreakQueue(network& _network, double _t) -> std::size_t;

private:
  auto breakData(const network& _network) -> std::tuple<double, size_t>;
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

private:
  auto readSavePoints(const toml::value& _config)
      -> quasiStaticStrainDouble::networkSavePoints;
};

}  // namespace protocols

}  // namespace networkV4