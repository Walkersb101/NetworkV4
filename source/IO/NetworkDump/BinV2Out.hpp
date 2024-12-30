#include <filesystem>
#include <fstream>
#include <optional>
#include <sstream>
#include <string>

#include <range/v3/view/iota.hpp>
#include <range/v3/view/zip.hpp>

#include "Core/Network.hpp"
#include "IO/BaseIO.hpp"
#include "IO/Compression/Compression.hpp"
#include "IO/NetworkDump/NetworkOut.hpp"

namespace IO
{
namespace networkDumps
{

class networkOutBinV2
    : public networkDump
    , public folderIO
    , public IO::Compression::compression
{
public:
  networkOutBinV2(const std::filesystem::path& _path,
                  IO::Compression::compressionTypes _compression,
                  const std::optional<size_t> _level = std::nullopt)
      : networkDump()
      , folderIO(_path)
      , IO::Compression::compression(_compression, _level) {};
  ~networkOutBinV2() {};

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override
  {
    auto fileName =
        _type + IO::Compression::compressionExts.lookup(m_type) + ".v2.bin";
    const auto path = m_dirPath / fileName;

    std::stringstream ss;

    const auto& nodes = _net.getNodes();
    const auto& bondMap = _net.getBonds().getMap();
    const auto& box = _net.getBox();

    append(ss, nodes.size());  // Append the size of the nodes
    append(ss, bondMap.size());  // Append the size of the bonds

    append(ss, box.getDomain());
    append(ss, box.shearStrain());

    for (const Utils::vec2d& pos : nodes.gatherPositions()) {
      append(ss, pos);
    }

    auto matrixTag = _net.getTags().get("matrix");

    for (const auto& [i, info, type, breakType] :
         ranges::views::zip(ranges::views::iota(0ul, bondMap.size()),
                            bondMap.getGlobalInfo(),
                            bondMap.getTypes(),
                            bondMap.getBreaks()))
    {
      if (!(std::holds_alternative<networkV4::Forces::VirtualBond>(type)
            || std::holds_alternative<networkV4::Forces::HarmonicBond>(
                type)))  // only two types of bonds are supported
      {
        throw("networkOutBinV2: unsupported bond type");
      }

      if (std::holds_alternative<networkV4::Forces::HarmonicBond>(type)
          && !std::get<networkV4::Forces::HarmonicBond>(type).normalized())
      {
        throw("networkOutBinV2: harmonic bond must be normalized");
      }

      append(ss, info.src);
      append(ss, info.dst);
      bool connected =
          std::holds_alternative<networkV4::Forces::HarmonicBond>(type);
      append(ss, connected);
      append(ss, bondMap.hasTag(i, matrixTag));
      if (connected) {
        append(ss, std::get<networkV4::Forces::HarmonicBond>(type).r0());
        append(ss, std::get<networkV4::Forces::HarmonicBond>(type).k());
        std::visit([&ss, this](const auto& _breakType) { append(ss, _breakType); },
                   breakType);
      } else {
        append(ss, 0.0);
        append(ss, 0.0);
        append(ss, 0.0);
      }
    }
    std::ofstream out(path, std::ios::binary);
    out << compress(ss.str());
    out.close();
  }

private:
  template<typename T>
  void append(std::stringstream& _ss, const T& _data) const
  {
    _ss.write(reinterpret_cast<const char*>(&_data), sizeof(T));
  }
};
}  // namespace networkDumps
}  // namespace IO
