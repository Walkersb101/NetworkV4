#pragma once

#include <filesystem>
#include <memory>

#include <toml.hpp>

#include "Core/Network.hpp"
#include "Misc/Enums.hpp"
#include "IO/NetworkOut.hpp"
#include "IO/DataOut.hpp"
//#include "IO/hdf5Out.hpp"
#include "Protocols/Protocol.hpp"
//#include "Protocols/StepStrain.hpp"
//#include "Protocols/Quasistatic.hpp"
//#include "Protocols/Propogator.hpp"

namespace networkV4
{
class Simulation
{
public:
  Simulation(const std::filesystem::path& _path);
  ~Simulation();


public:
  void run();

/*
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

*/
private:
  toml::value m_config;
  
  //networkV4::network m_network;
  std::unique_ptr<networkV4::protocol> m_protocol;
  std::unique_ptr<IO::timeSeries::dataOut> m_dataOut;
  std::unique_ptr<IO::networkDumps::networkOut> m_networkOut;
};
}  // namespace networkV4