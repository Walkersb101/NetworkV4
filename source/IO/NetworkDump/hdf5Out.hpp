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
    // initObservations(file, _net);
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

    /* TODO: Simulation Dependent
    // write observations
    auto sacrificialObs = file.getGroup("observables/sacrificialConnected");
    saveConnected(
        sacrificialObs, _net, _step, _time, networkV4::bondType::sacrificial);

    auto matrixObs = file.getGroup("observables/matrixConnected");
    saveConnected(matrixObs, _net, _step, _time, networkV4::bondType::matrix);
    */

    // write type
    // auto typeObs = file.getGroup("observables/logType");
    // saveScalar<int>(typeObs.getDataSet("step"), _step);
    // saveScalar<int>(typeObs.getDataSet("time"), _time);
    // saveScalar<std::string>(typeObs.getDataSet("value"), _type);

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

    std::vector<size_t> connections;
    connections.reserve(_net.getBonds().size() * 2);
    for (const networkV4::bonded::BondInfo& bond : _net.getBonds().getBonds()) {
      connections.emplace_back(bond.src);
      connections.emplace_back(bond.dst);
    }
    bonds.reshapeMemSpace({connections.size()}).write(connections);

    auto params = createOrGetGroup(_file, "parameter");
    auto bondParameters = initMatrixDataset<double>(params, "Bonds", {3, 2}, true);
    std::vector<double> data;
    data.reserve(_net.getBonds().size() * 3);
    for (auto [bond, type, brk] :
         ranges::views::zip(_net.getBonds().getBonds(),
                            _net.getBonds().getTypes(),
                            _net.getBonds().getBreaks()))
    {
        if (std::holds_alternative<networkV4::bonded::Forces::HarmonicBond>(type)) {
          const auto& bondType = std::get<networkV4::bonded::Forces::HarmonicBond>(type);
          data.push_back(bondType.mu());
          data.push_back(bondType.r0());
          data.push_back(bondType.lambda());
        } else {
          data.push_back(0.0);
          data.push_back(0.0);
          data.push_back(0.0);
        }
    }
    bondParameters.reshapeMemSpace({data.size()}).write(data);
  }
  // void initBondType(HighFive::File& _file,
  //                   const std::string& _name,
  //                   const networkV4::network& _net,
  //                   const networkV4::bondType _type) const;
  // void initObservations(HighFive::File& _file,
  //                       const networkV4::network& _net) const;

private:
  // void saveConnected(HighFive::Group& _group,
  //                    const networkV4::network& _net,
  //                    const std::size_t _step,
  //                    const double _time,
  //                    const networkV4::bondType _type) const;

private:
  // void saveConnectivity2Dataset(HighFive::DataSet& _dataset,
  //                               const networkV4::bonds& _bonds,
  //                               const networkV4::bondType _type) const;

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