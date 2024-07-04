#pragma once

#include <filesystem>
#include <memory>

#include <toml.hpp>

#include "Network.hpp"
#include "Protocol.hpp"
#include "DataOut.hpp"

class tomlIn
{
public:
  tomlIn(const std::filesystem::path& _path);
  ~tomlIn();

  auto readNetwork() -> networkV4::network;
  auto readProblem() -> std::unique_ptr<networkV4::protocol>;
  auto readDataOut() -> std::unique_ptr<dataOut>;

private:
  toml::value m_config;
};