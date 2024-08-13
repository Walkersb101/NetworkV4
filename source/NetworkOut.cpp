#include <fstream>
#include <iostream>
#include <stdexcept>

#include "NetworkOut.hpp"

#include <bxzstr.hpp>

#include "Bonds.hpp"
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

networkOut::networkOut(const std::filesystem::path& _path,
                       bxz::Compression _compression)
    : m_path(_path)
    , m_compression(_compression)
{
  if (!std::filesystem::exists(m_path.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_path.parent_path().string());
  }
}

void networkOut::save(const networkV4::network& _net,
                      const std::string& _name) const
{
}

networkOut::~networkOut() {}

noNetworkOut::noNetworkOut() {}

noNetworkOut::~noNetworkOut() {}

void noNetworkOut::save(const networkV4::network& _net,
                        const std::string& _name) const
{
}

networkOutBinV2::networkOutBinV2(const std::filesystem::path& _path,
                                 bxz::Compression _compression)
    : networkOut(_path, _compression)
{
}

networkOutBinV2::~networkOutBinV2() {}

void networkOutBinV2::save(const networkV4::network& _net,
                           const std::string& _name) const
{
  const auto path = m_path
      / (_name + networkV4::enumString::compressionExt.at(m_compression)
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