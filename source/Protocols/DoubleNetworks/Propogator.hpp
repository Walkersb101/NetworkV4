#pragma once

#include <cstdint>

#include "Integration/Integrators/Adaptive.hpp"
#include "Integration/Minimizers/Fire2.hpp"
#include "Integration/Minimizers/AdaptiveHeunDecent.hpp"
#include "Misc/Config.hpp"
#include "Misc/Roots.hpp"
#include "Protocols/Protocol.hpp"
#include "Protocols/deform.hpp"
#include "Protocols/protocolReader.hpp"

namespace networkV4
{
namespace protocols
{

class propogatorDouble : public protocolBase
{
public:
  propogatorDouble() = delete;
  propogatorDouble(
      std::shared_ptr<deform::deformBase>& _deform,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
      std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
      std::shared_ptr<IO::networkDumps::networkDump>& _networkOut,
      const network& _network,
      const std::vector<double>& _strains = {},
      double _rootTol = config::rootMethods::targetTol,
      integration::AdaptiveParams _params = integration::AdaptiveParams(),
      const minimisation::minimiserParams& _minParams =
          minimisation::minimiserParams(),
      double _maxStep = config::protocols::maxStep);
  ~propogatorDouble() = default;

public:
  void run(network& _network) override;

private:
  void runLambda(network& _network);
  void runStrain(network& _network);

  void evalStrain(network& _network, double _targetStrain);
  void relax(network& _network);

  auto getMaxDataIndex(network& _network,
                       const Utils::Tags::tagFlags& _filter) -> size_t;
  void breakMostStrained(network& _network,
                         const Utils::Tags::tagFlags& _filter);
  auto breakData(const network& _network) -> std::tuple<double, size_t>;

  auto findSingleBreak(network& _network) -> bool;

private:
  auto genTimeData(const network& _network,
                   const std::string& _reason,
                   std::size_t _breakCount)
      -> std::vector<IO::timeSeries::writeableTypes>;
  auto genBondData(const network& _network, const breakInfo& _bond)
      -> std::vector<IO::timeSeries::writeableTypes>;

private:
  std::vector<double> m_strains = {};

  integration::AdaptiveParams m_params;
  minimisation::minimiserParams m_minParams;
  double m_rootTol;

  size_t m_strainCount = 0;
  double m_maxStep;
};

class propogatorDoubleReader : public protocolReader
{
public:
  propogatorDoubleReader() = default;
  ~propogatorDoubleReader() = default;

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
