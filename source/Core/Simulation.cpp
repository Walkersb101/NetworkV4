#include "Simulation.hpp"

// #include "Misc/EnumString.hpp"

//--------------------------------------------------------------------

networkV4::tomlLoad::tomlLoad(const std::filesystem::path& _path)
{
  checkload(_path);
  m_config = toml::parse(_path);

  readNeworkIN();
  readOutputConfig();
  readProtocol();
}

void networkV4::tomlLoad::readNeworkIN()
{
  if (!m_config.contains("LoadPath")) {
    throw std::runtime_error("No LoadPath found in config file");
  }
  std::filesystem::path networkPath =
      toml::find<std::string>(m_config, "LoadPath");

  const std::string version =
      toml::find_or<std::string>(m_config, "Version", "BinV2");

  if (version == "BinV1") {
    // TODO: Implement this
    throw std::runtime_error("Load version not implemented");
  } else if (version == "BinV2") {
    m_networkIn = std::make_unique<IO::NetworkIn::BinV2In>(networkPath);
  } else {
    throw std::runtime_error("Load version not implemented");
  }
}

void networkV4::tomlLoad::readProtocol()
{
  if (m_config.contains("QuasiStaticStrain")) {
    m_protocolReader =
        std::make_unique<protocols::quasiStaticStrainDoubleReader>();
  } else {
    throw std::runtime_error("No protocol found in config file");
  }
}

void networkV4::tomlLoad::readOutputConfig()
{
  if (!m_config.contains("Output")) {
    throw std::runtime_error("No Output found in config file");
  }
  auto outConfig = toml::find(m_config, "Output");

  if (outConfig.contains("OutputPath")) {
    m_timeSeriesPath = toml::find<std::string>(outConfig, "OutputPath");
    m_networkDumpPath = toml::find<std::string>(outConfig, "OutputPath");
  } else if (outConfig.contains("DataPath") && outConfig.contains("BondPath")) {
    m_timeSeriesPath = toml::find<std::string>(outConfig, "DataPath");
    m_networkDumpPath = toml::find<std::string>(outConfig, "BondPath");
  } else {
    throw std::runtime_error(
        "OutputPath value (or Both DataPath and BondPath) not found in config "
        "file");
  }
}

auto networkV4::tomlLoad::loadtimeSeries(const std::string& _name)
    -> std::shared_ptr<IO::timeSeries::timeSeriesOut>
{
  auto outConfig = toml::find(m_config, "Output");
  const std::string outString =
      toml::find_or<std::string>(outConfig, "OutputType", "CSV");

  if (outString == "CSV") {
    const std::filesystem::path path = m_timeSeriesPath / (_name + ".csv");
    return std::make_shared<IO::timeSeries::CSVOut>(path);
  } else {
    throw std::runtime_error("OutputType not implemented");
  }
}

auto networkV4::tomlLoad::loadNetworkOut(const std::string& _name,
                                         const network& _network)
    -> std::shared_ptr<IO::networkDumps::networkDump>
{
  auto outConfig = toml::find(m_config, "Output");
  const std::string outString =
      toml::find_or<std::string>(outConfig, "OutputType", "None");

  if (outString == "HDF5") {
    const std::string author =
        toml::find_or<std::string>(outConfig, "Author", "Unknown");
    const std::filesystem::path path = m_networkDumpPath / (_name + ".h5");
    return std::make_shared<IO::networkDumps::networkOutHDF5>(
        path, _network, author);

  } else if (outString == "BinV2") {
    const auto [compType, level] = loadCompression();
    const std::filesystem::path path = m_networkDumpPath / (_name + ".v2.bin");
    return std::make_shared<IO::networkDumps::networkOutBinV2>(
        path, compType, level);

  } else {
    return std::make_shared<IO::networkDumps::noNetworkOut>();
  }
}

void networkV4::tomlLoad::checkload(const std::filesystem::path& _path)
{
  if (!std::filesystem::exists(_path)) {
    throw std::runtime_error("File does not exist: " + _path.string());
  }
  std::ifstream file(_path, std::ios::in);
  if (!file.is_open()) {
    throw std::runtime_error("Cannot open file: " + _path.string());
  }
  file.close();
}

auto networkV4::tomlLoad::loadCompression()
    -> std::pair<IO::Compression::compressionTypes, std::optional<size_t>>
{
  auto outConfig = toml::find(m_config, "Output");

  if (!m_config.contains("Compression")) {
    return {IO::Compression::compressionTypes::None, std::nullopt};
  }

  const std::string compressionStr =
      toml::find<std::string>(outConfig, "Compression");

  IO::Compression::compressionTypes compType;
  if (compressionStr == "None") {
    compType = IO::Compression::compressionTypes::None;
  } else if (compressionStr == "Zstd") {
    compType = IO::Compression::compressionTypes::Zstd;
  } else {
    throw std::runtime_error("Unknown compression type");
  }

  std::optional<size_t> level = std::nullopt;
  if (outConfig.contains("Level")) {
    level = toml::find<size_t>(outConfig, "Level");
  }
  return {compType, level};
}

//----------------------------------------------------------------------------

networkV4::Simulation::Simulation(const std::filesystem::path& _path)
    : tomlLoad(_path)
    , m_network(m_networkIn->load())
{
  m_config = toml::parse(_path);
  // loadTypes();
  initNetwork();

  m_dataOut = loadtimeSeries("Data");
  m_bondsOut = loadtimeSeries("Breaks");
  m_networkOut = loadNetworkOut("NetworkDump", m_network);
  m_protocol = m_protocolReader->read(
      m_config, m_network, m_dataOut, m_bondsOut, m_networkOut);
}

networkV4::Simulation::~Simulation() {}

void networkV4::Simulation::run()
{
  m_protocol->run(m_network);
}

void networkV4::Simulation::initNetwork()
{
  auto& nodes = m_network.getNodes();
  auto& bonds = m_network.getBonds();
  const auto& box = m_network.getBox();

  partition::PartitionGenerator partGen;
  partGen.assignNodes(nodes.positions(), box);
  partGen.sortNodes(nodes);
  partGen.sortBonds(bonds, nodes);
  partGen.checkPasses(bonds);

auto test = bonds.gatherBonds();

#if defined(_OPENMP)
  OMP::threadPartitions = partGen.generatePartitions(nodes, bonds);
  OMP::passes = partGen.getPasses();
#endif

  m_network.computeForces();
}

/*
void networkV4::Simulation::readStepStrain()
{
  auto stepConfig = toml::find(m_config, "StepStrain");

  if (!stepConfig.contains("StrainType")) {
    throw std::runtime_error("StrainType not specified");
  }
  const std::string strainTypeStr =
      toml::find<std::string>(stepConfig, "StrainType");
  networkV4::StrainType strainType =
      enumString::str2StrainType.at(strainTypeStr);

  m_breakProtocol = readBreakProtocol();

  if (!stepConfig.contains("MaxStrain")) {
    throw std::runtime_error("MaxStrain not specified");
  }
  const double maxStrain = toml::find<double>(stepConfig, "MaxStrain");

  const double ESP = toml::find_or<double>(
      stepConfig, "ESP", config::integrators::adaptiveIntegrator::esp);
  const double timeScale = toml::find_or<double>(
      stepConfig, "TimeScale", config::protocols::stepStrain::timeScale);
  const double stressScale = toml::find_or<double>(
      stepConfig, "StressScale", config::protocols::stepStrain::stressScale);
  const double maxTime = toml::find<double>(stepConfig, "MaxTime");

  m_protocol = std::make_unique<networkV4::stepStrain>(maxStrain,
                                                       strainType,
                                                       m_breakProtocol,
                                                       maxTime,
                                                       ESP,
                                                       timeScale,
                                                       stressScale);
}

void networkV4::Simulation::readPropogator()
{
  auto propConfig = toml::find(m_config, "Propogator");

  if (!propConfig.contains("Strains")) {
    throw std::runtime_error("Strains not specified");
  }
  std::vector<double> strains =
      toml::find<std::vector<double>>(propConfig, "Strains");

  if (!propConfig.contains("StrainType")) {
    throw std::runtime_error("StrainType not specified");
  }
  const std::string strainTypeStr =
      toml::find<std::string>(propConfig, "StrainType");
  const networkV4::StrainType strainType =
      enumString::str2StrainType.at(strainTypeStr);

  const std::string breakTypeStr =
      toml::find_or<std::string>(propConfig, "BreakType", "Any");
  bondType bType = enumString::str2BondType.at(breakTypeStr);

  const double ESP = toml::find_or<double>(
      propConfig, "ESP", config::integrators::adaptiveIntegrator::esp);
  const double tol = toml::find_or<double>(
      propConfig, "Tol", config::integrators::miminizer::tol);
  const double maxStep = toml::find_or<double>(propConfig, "MaxStep", 0.0);

  const bool respectLambda =
      toml::find_or<bool>(propConfig, "RespectLambda", false);

  m_protocol = std::make_unique<networkV4::propogator>(
      strains, strainType, bType, ESP, tol, maxStep, respectLambda);
}
*/
