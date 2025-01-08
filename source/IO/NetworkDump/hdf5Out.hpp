#pragma once

#include <filesystem>
#include <span>
#include <vector>

#include <highfive/highfive.hpp>

#include "IO/BaseIO.hpp"
#include "Misc/Utils.hpp"
#include "NetworkOut.hpp"

namespace IO
{
namespace networkDumps
{

class networkOutHDF5
    : public networkDump
    , public fileIO
{
public:
  networkOutHDF5(const std::filesystem::path& _path,
                 const networkV4::network& _net,
                 const std::string& _author = "Unknown")
      : networkDump()
      , fileIO(_path)
  {
    HighFive::File file(_path,
                        HighFive::File::Overwrite);  // create or overwrite file

    writeHeader(file, _author);
    initParticles(file, _net);
    initConnectivity(file, _net);
    initObservations(file, _net);
  }
  ~networkOutHDF5() {};

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override
  {
    HighFive::File file(m_filePath, HighFive::File::ReadWrite);

    auto allParticles = file.getGroup("particles").getGroup("all");

    // write box
    auto edgesGroup = allParticles.getGroup("box").getGroup("edges");
    saveScalar(edgesGroup.getDataSet("step"), _step);
    saveScalar(edgesGroup.getDataSet("time"), _time);
    saveMatrix(edgesGroup.getDataSet("value"), _net.getBox().getBox());

    auto posView = Utils::spanView(
        _net.getNodes().positions());  // Todo: work out why span is not working
    std::vector<double> poses;
    poses.assign(posView.begin(), posView.end());

    // write positions
    auto positions = allParticles.getGroup("positions");
    saveScalar(positions.getDataSet("step"), _step);
    saveScalar(positions.getDataSet("time"), _time);
    saveMatrix(positions.getDataSet("value"), poses);

    // TODO: Simulation Dependent
    const auto types = _net.getBonds().gatherTypes();

    // write observations
    auto connectedGroup = file.getGroup("observables/connected");
    std::vector<bool> connected;
    connected.reserve(_net.getBonds().size());

    for (const auto& typ : types) {
      if (std::holds_alternative<networkV4::Forces::HarmonicBond>(typ)) {
        connected.push_back(true);
      } else {
        connected.push_back(false);
      }
    }
    saveScalar(connectedGroup.getDataSet("step"), _step);
    saveScalar(connectedGroup.getDataSet("time"), _time);
    saveVector(connectedGroup.getDataSet("value"), connected);

    // write type
    auto typeObs = file.getGroup("observables/logType");
    saveScalar<int>(typeObs.getDataSet("step"), _step);
    saveScalar<double>(typeObs.getDataSet("time"), _time);
    saveScalar<std::string>(typeObs.getDataSet("value"), _type);

    file.flush();
  }

private:
  void writeHeader(HighFive::File& _file, const std::string& _author) const
  {
    auto h5md = createOrGetGroup(_file, "h5md");
    // create version attribute
    h5md.createAttribute("version", std::array<int, 2> {1, 1});

    // create author group
    h5md.createGroup("author").createAttribute<std::string>("name", _author);

    // create creator group
    auto creator = h5md.createGroup("creator");
    creator.createAttribute<std::string>("name", "HighFive");
    creator.createAttribute<std::string>("version", "v3.0.0-beta1");
  }
  void initParticles(HighFive::File& _file, const networkV4::network& _net)
  {
    auto particles = createOrGetGroup(_file, "particles");
    auto allParticles = createOrGetGroup(particles, "all");

    auto box = createOrGetGroup(allParticles, "box");
    box.createAttribute("dimension", 2);
    box.createAttribute("boundary", std::array<bool, 2> {true, true});

    // create edges group
    auto edges = createOrGetGroup(box, "edges");
    initScalarDataset<int>(edges, "step");
    initScalarDataset<double>(edges, "time");
    initMatrixDataset<double>(edges, "value", {2, 2});

    // create positions group
    auto positions = createOrGetGroup(allParticles, "positions");
    initScalarDataset<int>(positions, "step");
    initScalarDataset<double>(positions, "time");
    initMatrixDataset<double>(positions, "value", {_net.getNodes().size(), 2});
  }
  void initConnectivity(HighFive::File& _file,
                        const networkV4::network& _net) const
  {
    HighFive::Group connectivity =
        createOrGetGroup(_file, "particles/all/connectivity");
    auto bonds = initMatrixDataset<int>(
        connectivity, "Bonds", {_net.getBonds().size(), 2}, true);
    bonds.createAttribute(
        "particles_group",
        HighFive::Reference(createOrGetGroup(_file, "particles"),
                            createOrGetGroup(_file, "particles/all")));

    const auto bondInfo = _net.getBonds().gatherBonds();
    const auto types = _net.getBonds().gatherTypes();
    const auto breaks = _net.getBonds().gatherBreaks();

    std::vector<size_t> connections;
    connections.reserve(_net.getBonds().size() * 2);
    for (const networkV4::bonded::BondInfo& bond : bondInfo) {
      connections.emplace_back(bond.src);
      connections.emplace_back(bond.dst);
    }
    bonds.reshapeMemSpace({connections.size()}).write(connections);

    auto params = createOrGetGroup(_file, "parameter");
    auto bondParameters = initMatrixDataset<double>(
        params, "Bonds", {_net.getBonds().size(), 3}, true);

    std::vector<double> data;
    data.reserve(_net.getBonds().size() * 3);
    for (auto [type, brk] : ranges::views::zip(types, breaks)) {
      if (std::holds_alternative<networkV4::Forces::HarmonicBond>(type)) {
        const auto& bondType = std::get<networkV4::Forces::HarmonicBond>(type);
        data.push_back(bondType.k());
        data.push_back(bondType.r0());
      } else {
        data.push_back(0.0);
        data.push_back(0.0);
      }
      if (std::holds_alternative<networkV4::BreakTypes::StrainBreak>(brk)) {
        const auto& breakType =
            std::get<networkV4::BreakTypes::StrainBreak>(brk);
        data.push_back(breakType.lambda());
      } else {
        data.push_back(0.0);
      }
    }
    bondParameters.reshapeMemSpace({data.size()}).write(data);

    auto map = _net.getTags().getVecs();

    auto tagGroup = createOrGetGroup(params, "tags");

    auto variable_stringtype = HighFive::VariableLengthStringType();
    tagGroup
        .createDataSet("values",
                       HighFive::DataSpace({_net.getTags().size()},
                                           {_net.getTags().size()}),
                       variable_stringtype)
        .write(map);

    auto bondTags = createDataSet<bool>(
        tagGroup,
        "bondTags",
        std::vector<std::size_t> {_net.getBonds().size(),
                                  _net.getTags().size()},
        std::vector<std::size_t> {_net.getBonds().size(),
                                  _net.getTags().size()},
        std::vector<hsize_t> {_net.getBonds().size(), _net.getTags().size()},
        true);
    bondTags.write(_net.getBonds().gatherTags());
  }
  void initObservations(HighFive::File _file,
                        const networkV4::network& _net) const
  {
    auto connected = createOrGetGroup(_file, "observables/connected");
    initScalarDataset<int>(connected, "step");
    initScalarDataset<double>(connected, "time");
    createDataSet<bool>(
        connected,
        "value",
        std::vector<std::size_t> {0, _net.getBonds().size()},
        std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED,
                                  _net.getBonds().size()},
        std::vector<hsize_t> {1, _net.getBonds().size()},
        true);

    auto logType = createOrGetGroup(_file, "observables/logType");
    initScalarDataset<int>(logType, "step");
    initScalarDataset<double>(logType, "time");
    createDataSet<std::string>(
        logType,
        "value",
        std::vector<std::size_t> {0},
        std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
        std::vector<hsize_t> {1},
        false);
  }

private:
  template<typename T>
  HighFive::Group createOrGetGroup(T& _file, const std::string& _name) const
  {
    if (_file.exist(_name)) {
      return _file.getGroup(_name);
    } else {
      return _file.createGroup(_name);
    }
  }

  template<typename T>
  void saveScalar(HighFive::DataSet _dataset, const T& _data) const
  {
    const auto dataSetDims = _dataset.getSpace().getDimensions();
    _dataset.resize({dataSetDims[0] + 1});
    _dataset.select({dataSetDims[0]}, {1}).write(_data);
  }

  template<typename T>
  void saveVector(HighFive::DataSet _dataset, const std::vector<T>& _data) const
  {
    const auto dataSetDims = _dataset.getSpace().getDimensions();
    _dataset.resize({dataSetDims[0] + 1, dataSetDims[1]});
    _dataset.select({dataSetDims[0], 0}, {1, dataSetDims[1]}).squeezeMemSpace({0}).write(_data);
  }

  template<typename T>
  void saveMatrix(HighFive::DataSet _dataset, const std::vector<T>& _data) const
  {
    const auto dataSetDims = _dataset.getSpace().getDimensions();
    _dataset.resize({dataSetDims[0] + 1, dataSetDims[1], dataSetDims[2]});
    _dataset.select({dataSetDims[0], 0, 0}, {1, dataSetDims[1], dataSetDims[2]})
        .squeezeMemSpace({0})
        .reshapeMemSpace({_data.size()})
        .write(_data);
  }

  template<typename T>
  void saveMatrix(HighFive::DataSet _dataset,
                  const std::vector<std::vector<T>>& _data) const
  {
    const auto dataSetDims = _dataset.getSpace().getDimensions();

    _dataset.resize({dataSetDims[0] + 1, dataSetDims[1], dataSetDims[2]});
    _dataset.select({dataSetDims[0], 0, 0}, {1, dataSetDims[1], dataSetDims[2]})
        .squeezeMemSpace({0})
        .write(_data);
  }

  template<typename T>
  auto createDataSet(HighFive::Group& _group,
                     const std::string& _name,
                     std::vector<std::size_t> _dims,
                     std::vector<std::size_t> _maxDims,
                     std::vector<hsize_t> _chunkSize,
                     bool _deflate) const -> HighFive::DataSet
  {
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(_chunkSize));
    if (_deflate)
      props.add(HighFive::Deflate(9));
    props.add(HighFive::Shuffle());
    return _group.createDataSet<T>(
        _name, HighFive::DataSpace(_dims, _maxDims), props);
  }

  template<typename T>
  auto initScalarDataset(HighFive::Group& _group,
                         const std::string& _name) const -> HighFive::DataSet
  {
    return createDataSet<T>(
        _group,
        _name,
        std::vector<std::size_t> {0},
        std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
        std::vector<hsize_t> {1},
        false);
  }

  template<typename T>
  auto initMatrixDataset(HighFive::Group& _group,
                         const std::string& _name,
                         const std::pair<size_t, size_t>& _size,
                         bool _timeIndependent = false) const
      -> HighFive::DataSet
  {
    if (_timeIndependent) {
      return createDataSet<T>(
          _group,
          _name,
          std::vector<std::size_t> {_size.first, _size.second},
          std::vector<std::size_t> {_size.first, _size.second},
          std::vector<hsize_t> {_size.first, _size.second},
          false);
    }
    return createDataSet<T>(
        _group,
        _name,
        std::vector<std::size_t> {0, _size.first, _size.second},
        std::vector<std::size_t> {
            HighFive::DataSpace::UNLIMITED, _size.first, _size.second},
        std::vector<hsize_t> {1, _size.first, _size.second},
        true);
  }
};

}  // namespace networkDumps
}  // namespace IO