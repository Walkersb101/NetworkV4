#pragma once

#include <filesystem>

#include "Core/Bonds.hpp"
#include "Core/Nodes.hpp"
#include "Core/box.hpp"
#include "IO/BaseIO.hpp"
#include "NetworkIn.hpp"

namespace IO
{
namespace NetworkIn
{

class BinV2In
    : public fileIO
    , public networkIn
{
public:
  BinV2In(const std::filesystem::path& _filePath)
      : fileIO(_filePath)
      , networkIn()
  {
  }
  ~BinV2In() override = default;

public:
  networkV4::network load() override
  {
    std::ifstream file(m_filePath, std::ios::in | std::ios::binary);
    if (!file) {
      throw std::runtime_error("File " + m_filePath.string()
                               + " cannot be opened.");
    }

    const size_t N = read<size_t>(file);
    const size_t B = read<size_t>(file);

    const Utils::Math::vec2d domain = read<Utils::Math::vec2d>(file);
    const double shearStrain = read<double>(file);
    const networkV4::box domainBox(domain, shearStrain * domain[1]);

    networkV4::network network(domainBox, N, B);

    networkV4::nodes& nodes = network.getNodes();
    for (size_t i = 0; i < N; ++i) {
      const Utils::Math::vec2d pos = read<Utils::Math::vec2d>(file);
      nodes.addNode(pos);
    }

    auto& tags = network.getTags();
    auto matrixTag = tags.add("matrix");
    auto sacrificationTag = tags.add("sacrificial");

    networkV4::bonded::bonds& bonds = network.getBonds();
    for (size_t i = 0; i < B; ++i) {
      const size_t index1 = read<size_t>(file);
      const size_t index2 = read<size_t>(file);
      const bool connected = read<bool>(file);
      const bool matrix = read<bool>(file);
      const double naturalLength = read<double>(file);
      const double constant = read<double>(file);
      const double lambda = read<double>(file);

      if (index1 >= N || index2 >= N) {
        throw std::runtime_error("Invalid bond indices: "
                                 + std::to_string(index1) + ", "
                                 + std::to_string(index2));
      }

      networkV4::bonded::bondTypes type;
      if (connected) {
        type = networkV4::Forces::HarmonicBond(constant, naturalLength, true);
      } else {
        type = networkV4::Forces::VirtualBond();
      }

      networkV4::bonded::breakTypes breakType;
      if (connected) {
        breakType = networkV4::BreakTypes::StrainBreak(lambda, naturalLength);
      } else {
        breakType = networkV4::BreakTypes::None();
      }

      size_t index = bonds.addBond(index1, index2, type, breakType, matrix ? matrixTag : sacrificationTag);
    }
    network.getStresses().init(matrixTag);
    network.getStresses().init(sacrificationTag);

    return network;
  };

  template<typename T>
  void read(std::ifstream& _file, T& _data)
  {
    _file.read(reinterpret_cast<char*>(&_data), sizeof(T));
  }

  template<typename T>
  T read(std::ifstream& _file)
  {
    T data;
    _file.read(reinterpret_cast<char*>(&data), sizeof(T));
    return data;
  }
};

}  // namespace NetworkIn
}  // namespace IO