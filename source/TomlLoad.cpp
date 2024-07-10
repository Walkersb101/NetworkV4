#include <fstream>
#include <iostream>
#include <stdexcept>

#include "TomlLoad.hpp"

#include "Config.hpp"
#include "DataOut.hpp"
#include "Enums.hpp"

void checkCanOpen(const std::filesystem::path& _filename)
{
  if (!std::filesystem::exists(_filename)) {
    throw std::runtime_error("File does not exist: " + _filename.string());
  }
  std::ifstream file(_filename, std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _filename.string());
  }
}

tomlIn::tomlIn(const std::filesystem::path& _path)
{
  checkCanOpen(_path);
  m_config = toml::parse(_path);
}

tomlIn::~tomlIn() {}

auto getVersion(const std::string& _version) -> networkV4::loadVersion
{
  if (_version == "binV1") {
    return networkV4::loadVersion::binV1;
  } else if (_version == "binV2") {
    return networkV4::loadVersion::binV2;
  } else {
    throw std::runtime_error("Unknown version: " + _version);
  }
}

auto tomlIn::readNetwork() -> networkV4::network
{
  if (!m_config.contains("LoadPath")) {
    throw std::runtime_error("No LoadPath found in config file");
  }
  std::filesystem::path netPath = toml::find<std::string>(m_config, "LoadPath");
  checkCanOpen(netPath);

  std::string version =
      toml::find_or<std::string>(m_config, "Version", "binV2");

  networkV4::network net;
  switch (getVersion(version)) {
    case networkV4::loadVersion::binV1: {
      std::ifstream file(netPath, std::ios::binary | std::ios::in);
      if (!m_config.contains("Lambda")) {
        throw std::runtime_error(
            "lambda (break threshold) needs to be specified when loading a "
            "network with version binV1");
      }
      double lambda = toml::find<double>(m_config, "Lambda");
      net.loadFromBinV1(file, lambda);
      break;
    }
    case networkV4::loadVersion::binV2: {
      std::ifstream file(netPath, std::ios::binary | std::ios::in);
      net.loadFromBinV2(file);
      break;
    }
    default:
      throw std::runtime_error("load version not implemented");
  }
  return net;
}

auto getProtocol(const toml::value& _config) -> networkV4::protocolType
{
  if (_config.contains("QuasiStaticStrain")) {
    return networkV4::protocolType::QuasisaticStrain;
  } else {
    throw std::runtime_error("No Protocol found in config file");
  }
}

auto getBreakType(const std::string& _str) -> networkV4::BreakType
{
  if (_str == "Strain") {
    return networkV4::BreakType::Strain;
  } else {
    return networkV4::BreakType::None;
  }
}

auto getStrainType(const std::string& _str) -> networkV4::StrainType
{
  if (_str == "Shear") {
    return networkV4::StrainType::Shear;
  } else if (_str == "Elongation") {
    return networkV4::StrainType::Elongation;
  } else {
    throw std::runtime_error("Unknown StrainType: " + _str);
  }
}

auto tomlIn::readProblem() -> std::unique_ptr<networkV4::protocol>
{
  networkV4::protocolType protoType = getProtocol(m_config);

  switch (protoType) {
    case networkV4::protocolType::QuasisaticStrain: {
      auto quasiConfig = toml::find(m_config, "QuasiStaticStrain");
      if (!quasiConfig.contains("MaxStrain")) {
        throw std::runtime_error("MaxStrain not specified");
      }
      double maxStrain = toml::find<double>(quasiConfig, "MaxStrain");
      if (!quasiConfig.contains("StrainType")) {
        throw std::runtime_error("StrainType not specified");
      }
      networkV4::StrainType strainType =
          getStrainType(toml::find<std::string>(quasiConfig, "StrainType"));
      networkV4::BreakType breakType = getBreakType(
          toml::find_or<std::string>(quasiConfig, "BreakType", "None"));
      double ESP = toml::find_or<double>(
          quasiConfig, "ESP", config::adaptiveIntergrator::esp);
      double tol =
          toml::find_or<double>(quasiConfig, "Tol", config::ITPMethod::tol);

      return std::make_unique<networkV4::quasiStaticStrain>(
          maxStrain, strainType, breakType, ESP, tol);
    }
    default:
      throw std::runtime_error("Protocol not implemented");
  }
}

auto getOutput(const std::string& _str) -> networkV4::dataOutType
{
  if (_str == "CSV") {
    return networkV4::dataOutType::CSV;
  } else {
    throw std::runtime_error("No OutputType found in config file");
  }
}

auto tomlIn::readDataOut() -> std::unique_ptr<dataOut>
{
  if (!m_config.contains("Output")) {
    throw std::runtime_error("No Output found in config file");
  }
  auto outConfig = toml::find(m_config, "Output");

  auto type = toml::find_or<std::string>(
      outConfig, "OutputType", config::IO::outputType);

  std::filesystem::path dataPath;
  std::filesystem::path bondPath;
  if (outConfig.contains("OutputPath")) {
    dataPath = toml::find<std::string>(outConfig, "OutputPath");
    bondPath = toml::find<std::string>(outConfig, "OutputPath");
  } else if (outConfig.contains("DataPath") && outConfig.contains("BondPath")) {
    dataPath = toml::find<std::string>(outConfig, "DataPath");
    bondPath = toml::find<std::string>(outConfig, "BondPath");
  } else {
    throw std::runtime_error(
        "No OutputPath or DataPath and BondPath found in "
        "config file");
  }

  std::string timeDataName = toml::find_or<std::string>(
      outConfig, "TimeDataName", config::IO::timeDataName);
  std::string bondDataName = toml::find_or<std::string>(
      outConfig, "BondDataName", config::IO::bondDataName);

  switch (getOutput(type)) {
    case networkV4::dataOutType::CSV: {
      return std::make_unique<CSVOut>(
          dataPath, bondPath, timeDataName, bondDataName);
    }
    default: {
      throw std::runtime_error("OutputType not implemented");
    }
  }
}

auto getNetworkOut(const std::string& _str) -> networkV4::networkOutType
{
  if (_str == "BinV2") {
    return networkV4::networkOutType::BinV2;
  } else {
    return networkV4::networkOutType::None;
  }
}

auto tomlIn::readNetworkOut() -> std::unique_ptr<networkOut>
{
  if (!m_config.contains("Output")) {
    throw std::runtime_error("No Output found in config file");
  }
  auto outConfig = toml::find(m_config, "Output");

  auto type = toml::find_or<std::string>(
      outConfig, "NetworkType", "None");

  std::filesystem::path path;
  if (outConfig.contains("NetworkPath")) {
    path = toml::find<std::string>(outConfig, "NetworkPath");
  } else if (outConfig.contains("OutputPath")) {
    path = toml::find<std::string>(outConfig, "OutputPath");
  } else {
    throw std::runtime_error(
        "No NetworkPath or OutputPath found in config file");
  }

  switch (getNetworkOut(type)) {
    case networkV4::networkOutType::BinV2: {
      return std::make_unique<networkOutBinV2>(path);
    }
    case networkV4::networkOutType::None: {
      return std::make_unique<noNetworkOut>();
    }
    default: {
      throw std::runtime_error("OutputType not implemented");
    }
  }
}