#include "Simulation.hpp"

// #include "Misc/EnumString.hpp"

networkV4::Simulation::Simulation(const std::filesystem::path& _path)
{
  IO::NetworkIn::BinV2In binv2(
      "/home/sam/Documents/Code/NetworkV4/test/DoubleTest.bin");
  network net = binv2.load();

  IO::networkDumps::networkOutHDF5 hdf5Out(
      "/home/sam/Documents/Code/NetworkV4/test/DoubleTest.h5", net);

  net.computeForces();
  net.setBox(box({25, 23 * 23 / 25}, 0.0));
  net.computeForces();

  auto stopfn = [](double t, double dt, network& net)
  {
    if (t > 1) {
      return true;
    }
    return false;
  };

  integration::Run<integration::AdaptiveOverdampedEulerHeun> run {
      integration::AdaptiveOverdampedEulerHeun()};

  //run.run(net, 10000, stopfn);

  hdf5Out.save(net, 0, 0.0, "test");

  partition::PartitionGenerator partGen;
  partGen.assignNodes(net.getNodes().positions(), net.getBox());
  partGen.sortNodes(net.getNodes());
  partGen.sortBonds(net.getBonds(), net.getNodes());
  auto parts = partGen.generatePartitions(net.getNodes(), net.getBonds());

  net.computeForces();

  net.computeForces(parts, partGen.getPasses());

  // tools::checkCanOpen(_path);
  // m_config = toml::parse(_path);
  //
  // loadTypes();
  // loadNetwork();
  // readDataOut();
  // readNetworkOut();
  // readRandom();
  //
  // readProtocol();
}

networkV4::Simulation::~Simulation() {}

void networkV4::Simulation::run()
{
  // m_protocol->initIO(m_network, m_dataOut, m_networkOut);
  // m_protocol->run(m_network);
}

/*
void networkV4::Simulation::loadTypes()
{
  if (!m_config.contains("LoadPath")) {
    throw std::runtime_error("No LoadPath found in config file");
  }
  m_networkPath = toml::find<std::string>(m_config, "LoadPath");
  tools::checkCanOpen(m_networkPath);

  std::string version =
      toml::find_or<std::string>(m_config, "Version", "BinV2");
  m_loadVersion = enumString::str2LoadVersion.at(version);

  if (!m_config.contains("Output")) {
    throw std::runtime_error("No Output found in config file");
  }
  auto outConfig = toml::find(m_config, "Output");
  const std::string outString = toml::find_or<std::string>(
      outConfig, "OutputType", config::IO::outputType);
  m_dataOutType = enumString::str2DataOutType.at(outString);
  const std::string typeStr =
      toml::find_or<std::string>(outConfig, "NetworkType", "None");
  m_networkOutType = enumString::str2NetworkOutType.at(typeStr);

  m_protocolType = getProtocol();
  m_breakType = networkV4::BreakType::None;
}

void networkV4::Simulation::loadNetwork()
{
  std::ifstream file(m_networkPath, std::ios::binary | std::ios::in);

  switch (m_loadVersion) {
    case networkV4::loadVersion::BinV1: {
      if (!m_config.contains("Lambda")) {
        throw std::runtime_error(
            "lambda (break threshold) needs to be specified when loading a "
            "network with version BinV1");
      }
      const double lambda = toml::find<double>(m_config, "Lambda");
      m_network.loadFromBinV1(file, lambda);
      break;
    }
    case networkV4::loadVersion::BinV2: {
      m_network.loadFromBinV2(file);
      break;
    }
    default:
      throw std::runtime_error("load version not implemented");
  }
  m_network.computeForces();
}

void networkV4::Simulation::readDataOut()
{
  auto outConfig = toml::find(m_config, "Output");

  std::filesystem::path dataPath;
  std::filesystem::path statesPath;
  if (outConfig.contains("OutputPath")) {
    dataPath = toml::find<std::string>(outConfig, "OutputPath");
    statesPath = toml::find<std::string>(outConfig, "OutputPath");
  } else if (outConfig.contains("DataPath") && outConfig.contains("BondPath"))
{ dataPath = toml::find<std::string>(outConfig, "DataPath"); statesPath =
toml::find<std::string>(outConfig, "BondPath"); } else { throw
std::runtime_error( "OutputPath value (or Both DataPath and BondPath found in
" "config file");
  }

  const std::string timeDataName = toml::find_or<std::string>(
      outConfig, "TimeDataName", config::IO::timeDataName);
  const std::string stateDataName = toml::find_or<std::string>(
      outConfig, "BondDataName", config::IO::bondDataName);

  switch (m_dataOutType) {
    case networkV4::dataOutType::CSV: {
      m_dataOut = std::make_unique<CSVOut>(
          dataPath, statesPath, timeDataName, stateDataName);
      break;
    }
    default: {
      throw std::runtime_error("OutputType not implemented");
    }
  }
}

void networkV4::Simulation::readNetworkOut()
{
  auto outConfig = toml::find(m_config, "Output");

  const std::string compressionStr =
      toml::find_or<std::string>(outConfig, "Compression", "PlainText");
  const bxz::Compression compression =
      enumString::str2Compression.at(compressionStr);

  std::filesystem::path path;
  if (outConfig.contains("NetworkPath")) {
    path = toml::find<std::string>(outConfig, "NetworkPath");
  } else if (outConfig.contains("OutputPath")) {
    path = toml::find<std::string>(outConfig, "OutputPath");
  } else {
    throw std::runtime_error(
        "No NetworkPath or OutputPath found in config file");
  }

  switch (m_networkOutType) {
    case networkV4::networkOutType::HDF5: {
      m_networkOut = std::make_unique<networkOutHDF5>(path, m_network);
      break;
    }
    case networkV4::networkOutType::BinV2: {
      m_networkOut = std::make_unique<networkOutBinV2>(path, compression);
      break;
    }
    case networkV4::networkOutType::None: {
      m_networkOut = std::make_unique<noNetworkOut>();
      break;
    }
    default: {
      throw std::runtime_error("OutputType not implemented");
    }
  }
}

void networkV4::Simulation::readRandom()
{
  if (!m_config.contains("Random")) {
    return;
  }
  auto randomConfig = toml::find(m_config, "Random");
  if (randomConfig.contains("Seed")) {
    m_seed = toml::find<unsigned long>(randomConfig, "Seed");
  } else {
    m_seed =
        (unsigned long)std::chrono::duration_cast<std::chrono::milliseconds>(
            std::chrono::system_clock::now().time_since_epoch())
            .count();
  }
}

auto networkV4::Simulation::loadBreakType() -> networkV4::BreakType
{
  std::string breakTypeStr;
  toml::value config;
  if (m_config.contains("BreakType")) {
    config = toml::find(m_config, "BreakType");
    breakTypeStr = toml::find_or<std::string>(config, "Type", "None");
  } else {
    toml::value protoConfig =
        toml::find(m_config, enumString::protocolType2Str.at(m_protocolType));
    if (protoConfig.contains("BreakType")) {
      breakTypeStr =
          toml::find_or<std::string>(protoConfig, "BreakType", "None");
    } else {
      throw std::runtime_error("No BreakType found in config file");
    }
  }
  return enumString::str2BreakType.at(breakTypeStr);
}

auto networkV4::Simulation::readBreakProtocol()
    -> std::unique_ptr<networkV4::BreakTypes>
{
  toml::value config =
      toml::find(m_config, enumString::protocolType2Str.at(m_protocolType));
  std::string breakTypeStr =
      toml::find_or<std::string>(config, "BreakType", "None");
  m_breakType = enumString::str2BreakType.at(breakTypeStr);
  switch (m_breakType) {
    case networkV4::BreakType::None: {
      return std::make_unique<networkV4::NoBreaks>();
    }
    case networkV4::BreakType::Strain: {
      std::vector<double> lambda = m_network.getBonds().getLambdaVector();
      return std::make_unique<networkV4::StrainBreak>();
    }
    case networkV4::BreakType::Energy: {
      std::vector<double> lambda;
      if (config.contains("E")) {
        double E = toml::find<double>(config, "E");
        lambda = std::vector<double>(m_network.getBonds().size(), E);
        m_network.getBonds().setLambdaVector(lambda);
      }
      return std::make_unique<networkV4::EnergyBreak>();
    }
    case networkV4::BreakType::SGR: {
      std::vector<double> lambda;
      if (config.contains("E")) {
        double E = toml::find<double>(config, "E");
        lambda = std::vector<double>(m_network.getBonds().size(), E);
        m_network.getBonds().setLambdaVector(lambda);
      }
      if (!config.contains("Tau0")) {
        throw std::runtime_error("Tao0 not specified");
      }
      double tau0 = toml::find_or<double>(config, "Tau0", 1.0);
      if (!config.contains("T")) {
        throw std::runtime_error("T not specified");
      }
      double T = toml::find<double>(config, "T");
      return std::make_unique<networkV4::SGRBreak>(tau0, T, m_seed);
    }
    default: {
      throw std::runtime_error("BreakType not implemented");
    }
  }
}

auto networkV4::Simulation::getProtocol() -> networkV4::protocolType
{
  for (const auto& [key, value] : enumString::str2ProtocolType) {
    if (m_config.contains(key)) {
      return value;
    }
  }
  throw std::runtime_error("No Protocol found in config file");
}

void networkV4::Simulation::readProtocol()
{
  switch (m_protocolType) {
    case networkV4::protocolType::QuasisaticStrain: {
      readQuasistatic();
      break;
    }
    case networkV4::protocolType::StepStrain: {
      readStepStrain();
      break;
    }
    case networkV4::protocolType::Propogator: {
      readPropogator();
      break;
    }
    default:
      throw std::runtime_error("Protocol not implemented");
  }
}

void networkV4::Simulation::readQuasistatic()
{
  auto quasiConfig = toml::find(m_config, "QuasiStaticStrain");

  if (!quasiConfig.contains("MaxStrain")) {
    throw std::runtime_error("MaxStrain not specified");
  }
  const double maxStrain = toml::find<double>(quasiConfig, "MaxStrain");

  if (!quasiConfig.contains("StrainType")) {
    throw std::runtime_error("StrainType not specified");
  }
  const std::string strainTypeStr =
      toml::find<std::string>(quasiConfig, "StrainType");
  networkV4::StrainType strainType =
      enumString::str2StrainType.at(strainTypeStr);

  m_breakProtocol = readBreakProtocol();

  const double ESP = toml::find_or<double>(
      quasiConfig, "ESP", config::integrators::adaptiveIntegrator::esp);
  const double tol =
      toml::find_or<double>(quasiConfig, "Tol",
config::rootMethods::targetTol); const bool errorOnNotSingleBreak =
toml::find_or<bool>( quasiConfig, "ErrorOnNotSingleBreak",
      config::protocols::quasiStaticStrain::errorOnNotSingleBreak);
  const double maxStep =
      toml::find_or<double>(quasiConfig, "MaxStep", maxStrain);
  const bool saveBreaks = toml::find_or<bool>(quasiConfig, "SaveBreaks",
false);

  m_protocol =
      std::make_unique<networkV4::quasiStaticStrain>(maxStrain,
                                                     strainType,
                                                     m_breakProtocol,
                                                     ESP,
                                                     tol,
                                                     errorOnNotSingleBreak,
                                                     maxStep,
                                                     saveBreaks);
}

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

  const bool respectLambda = toml::find_or<bool>(
      propConfig, "RespectLambda", false);

  m_protocol = std::make_unique<networkV4::propogator>(
      strains, strainType, bType, ESP, tol, maxStep, respectLambda);
}
*/