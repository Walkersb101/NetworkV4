#include <fstream>
#include <iostream>
#include <stdexcept>

#include "NetworkOut.hpp"

#include <bxzstr.hpp>
#include <highfive/highfive.hpp>

#include "Bonds.hpp"
#include "Config.hpp"
#include "Network.hpp"
#include "Nodes.hpp"
#include "Tools.hpp"

template<>
bool write<bxz::Compression>(const std::stringstream& _ss,
                             const std::filesystem::path& _path,
                             bxz::Compression _compression)
{
  const auto mode = std::ios::out | std::ios::binary;
  if (_compression == bxz::Compression::plaintext) {
    return write(_ss, _path, std::ofstream(_path, mode));
  } else {
    return write(_ss, _path, bxz::ofstream(_path, mode, _compression));
  }
}

networkOut::networkOut() {}

networkOut::networkOut(const std::filesystem::path& _path)
    : m_path(_path)
{
  if (!std::filesystem::exists(m_path.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_path.parent_path().string());
  }
}

void networkOut::save(const networkV4::network& _net,
                      const std::size_t _step,
                      const double _time,
                      const std::string& _type) const
{
}

networkOut::~networkOut() {}

noNetworkOut::noNetworkOut() {}

noNetworkOut::~noNetworkOut() {}

void noNetworkOut::save(const networkV4::network& _net,
                        const std::size_t _step,
                        const double _time,
                        const std::string& _type) const
{
}

networkOutBinV2::networkOutBinV2(const std::filesystem::path& _path,
                                 bxz::Compression _compression)
    : networkOut(_path)
    , m_compression(_compression)
{
}

networkOutBinV2::~networkOutBinV2() {}

void networkOutBinV2::save(const networkV4::network& _net,
                           const std::size_t _step,
                           const double _time,
                           const std::string& _type) const
{
  const auto path = m_path
      / (_type + networkV4::enumString::compressionExt.at(m_compression)
         + ".v2.bin");

  std::stringstream ss;

  const networkV4::nodes& nodes = _net.getNodes();
  const networkV4::bonds& bonds = _net.getBonds();

  const std::size_t N = nodes.size();
  const std::size_t B = bonds.size();
  ss.write(reinterpret_cast<const char*>(&N), sizeof(std::size_t));
  ss.write(reinterpret_cast<const char*>(&B), sizeof(std::size_t));

  const vec2d& domain = _net.getDomain();
  const double shearStrain = _net.getShearStrain();
  ss.write(reinterpret_cast<const char*>(&domain), sizeof(vec2d));
  ss.write(reinterpret_cast<const char*>(&shearStrain), sizeof(double));

  for (const vec2d& node : nodes.positions()) {
    ss.write(reinterpret_cast<const char*>(&node), sizeof(vec2d));
  }

  for (const networkV4::bond& b : bonds) {
    const std::size_t index1 = b.src();
    const std::size_t index2 = b.dst();
    const bool connected = b.connected();
    const bool matrix = b.type() == networkV4::bondType::matrix
        || b.type() == networkV4::bondType::single;
    const double naturalLength = b.naturalLength();
    const double constant = b.mu();
    const double lambda = b.lambda();

    ss.write(reinterpret_cast<const char*>(&index1), sizeof(std::size_t));
    ss.write(reinterpret_cast<const char*>(&index2), sizeof(std::size_t));
    ss.write(reinterpret_cast<const char*>(&connected), sizeof(bool));
    ss.write(reinterpret_cast<const char*>(&matrix), sizeof(bool));
    ss.write(reinterpret_cast<const char*>(&naturalLength), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&constant), sizeof(double));
    ss.write(reinterpret_cast<const char*>(&lambda), sizeof(double));
  }
  write(ss, path, m_compression);
  ss.clear();
}

networkOutHDF5::networkOutHDF5(const std::filesystem::path& _path)
    : networkOut(_path)
{
  HighFive::DataSetCreateProps props;

  // create or overwrite file
  m_file = HighFive::File(m_path.string(), HighFive::File::Overwrite);

  // create h5md group
  auto h5md = m_file.createGroup("h5md");

  // create version attribute
  h5md.createAttribute("version", std::array<int, 2> {1, 1});

  // create author group
  // TODO: add author
  auto author = h5md.createGroup("author");
  author.createAttribute("name", "Network V4");

  // create creator group
  auto creator = h5md.createGroup("creator");
  creator.createAttribute("name", "HighFive");
  creator.createAttribute("version", "v3.0.0-beta1");

  // create positions group
  auto particles = h5md.createGroup("particles");
  auto allParticles = h5md.createGroup("all");

  // setup box
  auto box = allParticles.createGroup("box");
  box.createAttribute("dimension", 2);
  box.createAttribute("boundary", std::array<bool, 2> {true, true});
  auto edges = box.createGroup("edges");

  // Scalar Properties
  props = HighFive::DataSetCreateProps();
  props.add(HighFive::Chunking(std::vector<hsize_t> {1}));

  edges.createDataSet<int>(
      "step",
      HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}),
      props);
  edges.createDataSet<double>(
      "time",
      HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}),
      props);

  // Vector Properties
  props = HighFive::DataSetCreateProps();
  props.add(HighFive::Chunking(std::vector<hsize_t> {1, 2}));

  edges.createDataSet<double>(
      "value",
      HighFive::DataSpace({0, 2}, {HighFive::DataSpace::UNLIMITED, 2}),
      props);

  // create positions group
  auto positions = allParticles.createGroup("positions");

  props = HighFive::DataSetCreateProps();
  props.add(HighFive::Chunking(std::vector<hsize_t> {1}));

  positions.createDataSet<int>(
      "step",
      HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}),
      props);
  positions.createDataSet<double>(
      "time",
      HighFive::DataSpace({0}, {HighFive::DataSpace::UNLIMITED}),
      props);

  props = HighFive::DataSetCreateProps();
  props.add(
      HighFive::Chunking(std::vector<hsize_t> {config::hdf5::maxChunk / 2, 2}));
}

networkOutHDF5::~networkOutHDF5() {}

void networkOutHDF5::save(const networkV4::network& _net,
                           const std::size_t _step,
                           const double _time,
                           const std::string& _type) const
{
}