#pragma once

#include <filesystem>
#include <memory>

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Enums.hpp"
#include "IO/NetworkOut.hpp"
#include "IO/DataOut.hpp"
#include "IO/hdf5Out.hpp"
#include "Protocols/Protocol.hpp"
#include "Protocols/StepStrain.hpp"
#include "Protocols/Quasistatic.hpp"
#include "Protocols/Propogator.hpp"

namespace networkV4
{
class Simulation
{
public:
  Simulation(const std::filesystem::path& _path);
  ~Simulation();

public:
  void run();

private:
  void loadTypes();

  void loadNetwork();
  void readDataOut();
  void readNetworkOut();
  void readRandom();

  auto loadBreakType() -> networkV4::BreakType;
  auto readBreakProtocol() -> std::unique_ptr<BreakTypes>;

private:
  auto getProtocol() -> networkV4::protocolType;
  void readProtocol();
  void readQuasistatic();
  void readStepStrain();
  void readPropogator();

private:
  toml::value m_config;

  std::filesystem::path m_networkPath;

  networkV4::loadVersion m_loadVersion;
  networkV4::dataOutType m_dataOutType;

  networkV4::networkOutType m_networkOutType;
  networkV4::protocolType m_protocolType;
  
  networkV4::BreakType m_breakType;

  networkV4::network m_network;
  std::unique_ptr<networkV4::protocol> m_protocol;
  std::unique_ptr<dataOut> m_dataOut;
  std::unique_ptr<networkOut> m_networkOut;
  std::unique_ptr<BreakTypes> m_breakProtocol;

  unsigned long m_seed;
};
}  // namespace networkV4