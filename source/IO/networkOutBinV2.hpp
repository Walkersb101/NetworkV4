#include "IO/NetworkOut.hpp"

namespace IO
{
namespace networkDumps
{

template<typename T>
bool write(const std::stringstream& _ss,
           const std::filesystem::path& _path,
           T _file)
{
  if (!_file.is_open()) {
    return false;
  }
  _file << _ss.str() << std::flush;
  _file.close();
  return _file.good();
}

class networkOutBinV2 : public networkOut
{
public:
  networkOutBinV2(const std::filesystem::path& _path,
                  Compression _compression)
      : networkOut(_path)
      , m_compression(_compression)
  {
  }
  ~networkOutBinV2();

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override
  {
    const auto path = m_path
        / (_type + networkV4::enumString::compressionExt.at(m_compression)
           + ".v2.bin");

    std::stringstream ss;

    const networkV4::nodes& nodes = _net.getNodes();
    const networkV4::bonded::bonds& bonds = _net.getBonds();

    const std::size_t N = nodes.size();
    const std::size_t B = bonds.size();
    ss.write(reinterpret_cast<const char*>(&N), sizeof(std::size_t));
    ss.write(reinterpret_cast<const char*>(&B), sizeof(std::size_t));

    const Utils::vec2d& domain = _net.getBox().getDomain();
    const double shearStrain = _net.getBox().shearStrain();
    ss.write(reinterpret_cast<const char*>(&domain), sizeof(Utils::vec2d));
    ss.write(reinterpret_cast<const char*>(&shearStrain), sizeof(double));

    for (const Utils::vec2d& node : nodes.positions()) {
      ss.write(reinterpret_cast<const char*>(&node), sizeof(Utils::vec2d));
    }

    /*
    for (const networkV4::bonded::& b : bonds.getMap().getGlobalInfo())
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
    */ //TODO: Fix this
    write(ss, path, m_compression);
    ss.clear();
  }

private:
  bxz::Compression m_compression;
};
}  // namespace networkDumps
}  // namespace IO