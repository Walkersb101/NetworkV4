#pragma once

#include <toml.hpp>

#include "Integration/Minimizers/MinimiserBase.hpp"
#include "Protocols/Protocol.hpp"

namespace networkV4
{

namespace protocols
{
class protocolReader
{
public:
  protocolReader() = default;
  virtual ~protocolReader() = default;

public:
  virtual auto read(const toml::value& _config,
                    const network& _network,
                    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _dataOut,
                    std::shared_ptr<IO::timeSeries::timeSeriesOut>& _bondsOut,
                    std::shared_ptr<IO::networkDumps::networkDump>& _networkOut)
      -> std::shared_ptr<protocolBase> = 0;

protected:
  auto readDeform(const toml::value& _config)
      -> std::shared_ptr<deform::deformBase>
  {
    if (!_config.contains("StrainType"))
      throw std::runtime_error("StrainType not specified");
    const std::string deformType =
        toml::find<std::string>(_config, "StrainType");
    if (deformType == "Shear")
      return std::make_shared<deform::shear>();
    if (deformType == "Elongation")
      return std::make_shared<deform::elongationAreaY>();
    throw std::runtime_error("Deform type not implemented");
  }

  auto readAdaptive(const toml::value& _config) -> integration::AdaptiveParams
  {
    if (!_config.contains("Adaptive"))
      return {};
    auto config = toml::find(_config, "Adaptive");
    integration::AdaptiveParams params;
    if (config.contains("ESP"))
      params.espAbs = toml::find<double>(config, "ESPabs");
    if (config.contains("ESPrel"))
      params.espRel = toml::find<double>(config, "ESPrel");
    if (config.contains("qMin"))
      params.qMin = toml::find<double>(config, "qMin");
    if (config.contains("qMax"))
      params.qMax = toml::find<double>(config, "qMax");
    if (config.contains("dtMin"))
      params.dtMin = toml::find<double>(config, "dtMin");
    if (config.contains("dtMax"))
      params.dtMax = toml::find<double>(config, "dtMax");
    if (config.contains("maxInnerIter"))
      params.maxInnerIter = toml::find<size_t>(config, "maxInnerIter");

    return params;
  }

  auto readMinimiser(const toml::value& _config)
      -> minimisation::minimiserParams
  {
    if (!_config.contains("Minimiser"))
      return {};
    auto config = toml::find(_config, "Minimiser");
    minimisation::minimiserParams params;
    if (config.contains("Ftol"))
      params.Ftol = toml::find<double>(config, "Ftol");
    if (config.contains("Etol"))
      params.Etol = toml::find<double>(config, "Etol");
    if (config.contains("maxIter"))
      params.maxIter = toml::find<size_t>(config, "maxIter");

    return params;
  }
};

}  // namespace protocols
}  // namespace networkV4