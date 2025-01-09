#pragma once

#include <filesystem>
#include <memory>

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Core/OMP/OMP.hpp"
#include "IO/Input/Binv2.hpp"
#include "IO/NetworkDump/BinV2Out.hpp"
#include "IO/NetworkDump/NetworkOut.hpp"
#include "IO/NetworkDump/hdf5Out.hpp"
#include "IO/TimeSeries/CSVOut.hpp"
#include "IO/TimeSeries/DataOut.hpp"
#include "Integration/RunStyles/Run.hpp"
#include "Misc/Config.hpp"
#include "Protocols/Protocol.hpp"
#include "Protocols/protocolReader.hpp"
#include "Protocols/DoubleNetworks/Quasistatic.hpp"
#include "Protocols/DoubleNetworks/Propogator.hpp"

namespace networkV4
{

class tomlLoad
{
public:
  tomlLoad() = delete;
  tomlLoad(const std::filesystem::path& _path);
  ~tomlLoad() = default;

private:
  void readNeworkIN();
  void readProtocol();

protected:
  void readOutputConfig();
  auto loadtimeSeries(const std::string& _name)
      -> std::shared_ptr<IO::timeSeries::timeSeriesOut>;
  auto loadNetworkOut(const std::string& _name, const network& _network)
      -> std::shared_ptr<IO::networkDumps::networkDump>;

protected:
  toml::value m_config;

  std::unique_ptr<IO::NetworkIn::networkIn> m_networkIn;
  std::unique_ptr<protocols::protocolReader> m_protocolReader;
  std::filesystem::path m_timeSeriesPath;
  std::filesystem::path m_networkDumpPath;

private:
  void checkload(const std::filesystem::path& _path);
  auto loadCompression()
      -> std::pair<IO::Compression::compressionTypes, std::optional<size_t>>;
};

class Simulation : public tomlLoad
{
public:
  Simulation(const std::filesystem::path& _path);
  ~Simulation();

public:
  void run();

private:
  // void loadTypes();

  void initNetwork();
  // void readDataOut();
  // void readNetworkOut();
  // void readRandom();

  // auto loadBreakType() -> networkV4::BreakType;
  // auto readBreakProtocol() -> std::unique_ptr<BreakTypes>;

private:
  // auto getProtocol() -> networkV4::protocolType;
  // void readProtocol();
  // void readQuasistatic();
  // void readStepStrain();
  // void readPropogator();

private:
  toml::value m_config;

  networkV4::network m_network;
  std::shared_ptr<networkV4::protocols::protocolBase> m_protocol;
  std::shared_ptr<IO::timeSeries::timeSeriesOut> m_dataOut;
  std::shared_ptr<IO::timeSeries::timeSeriesOut> m_bondsOut;
  std::shared_ptr<IO::networkDumps::networkDump> m_networkOut;
};
}  // namespace networkV4