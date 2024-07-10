#include <fstream>
#include <iostream>
#include <stdexcept>

#include "NetworkOut.hpp"

#include "Bonds.hpp"
#include "Network.hpp"
#include "Nodes.hpp"

networkOut::networkOut() {}

networkOut::networkOut(const std::filesystem::path& _path)
    : m_path(_path)
{
  if (!std::filesystem::exists(m_path.parent_path())) {
    throw std::runtime_error("Directory does not exist: "
                             + m_path.parent_path().string());
  }
}

networkOut::~networkOut() {}

noNetworkOut::noNetworkOut() {}

noNetworkOut::~noNetworkOut() {}

void noNetworkOut::save(const networkV4::network& _net, const std::string& _name) const {}

networkOutBinV2::networkOutBinV2(const std::filesystem::path& _path)
    : networkOut(_path)
{
}

networkOutBinV2::~networkOutBinV2() {}

void networkOutBinV2::save(const networkV4::network& _net, const std::string& _name) const
{
  std::ofstream file(m_path / _name, std::ios_base::out | std::ios_base::binary);
  if (!file.is_open()) {
    std::runtime_error("Could not open file for writing");
  }

  const networkV4::nodes& nodes = _net.getNodes();
  const networkV4::bonds& bonds = _net.getBonds();

  const std::size_t N = nodes.size();
  const std::size_t B = bonds.size();
  file.write(reinterpret_cast<const char*>(&N), sizeof(std::size_t));
  file.write(reinterpret_cast<const char*>(&B), sizeof(std::size_t));

  const vec2d& domain = _net.getDomain();
  const double shearStrain = _net.getShearStrain();
  file.write(reinterpret_cast<const char*>(&domain), sizeof(vec2d));
  file.write(reinterpret_cast<const char*>(&shearStrain), sizeof(double));

  for (const vec2d& node : nodes.positions()) {
    file.write(reinterpret_cast<const char*>(&node), sizeof(vec2d));
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

    file.write(reinterpret_cast<const char*>(&index1), sizeof(std::size_t));
    file.write(reinterpret_cast<const char*>(&index2), sizeof(std::size_t));
    file.write(reinterpret_cast<const char*>(&connected), sizeof(bool));
    file.write(reinterpret_cast<const char*>(&matrix), sizeof(bool));
    file.write(reinterpret_cast<const char*>(&naturalLength), sizeof(double));
    file.write(reinterpret_cast<const char*>(&constant), sizeof(double));
    file.write(reinterpret_cast<const char*>(&lambda), sizeof(double));
  }
  file.close();
}