#pragma once

#include <filesystem>

#include <highfive/highfive.hpp>

#include "BaseIO.hpp"
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
    saveBox(edgesGroup, _net, _step, _time);

    // write positions
    auto positions = allParticles.getGroup("positions");
    savePositions(positions, _net, _step, _time);

    /* TODO: Simulation Dependent
    // write observations
    auto sacrificialObs = file.getGroup("observables/sacrificialConnected");
    saveConnected(
        sacrificialObs, _net, _step, _time, networkV4::bondType::sacrificial);

    auto matrixObs = file.getGroup("observables/matrixConnected");
    saveConnected(matrixObs, _net, _step, _time, networkV4::bondType::matrix);
    */

    // write type
    auto typeObs = file.getGroup("observables/logType");
    saveScalar<int>(typeObs.getDataSet("step"), _step);
    saveScalar<int>(typeObs.getDataSet("time"), _time);
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

    initBox(allParticles);

    // create positions group
    auto positions = createOrGetGroup(allParticles, "positions");
    size_t posChunkSize =
        findChunkSize(_net.getNodes().size(), config::hdf5::maxChunk / 2);

    createDataSet<int>(
        positions,
        "step",
        std::vector<std::size_t> {0},
        std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
        std::vector<hsize_t> {1},
        false);
    createDataSet<double>(
        positions,
        "time",
        std::vector<std::size_t> {0},
        std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
        std::vector<hsize_t> {1},
        false);

    createDataSet<double>(
        positions,
        "value",
        std::vector<std::size_t> {0, _net.getNodes().size(), 2},
        std::vector<std::size_t> {
            HighFive::DataSpace::UNLIMITED, _net.getNodes().size(), 2},
        std::vector<hsize_t> {1, posChunkSize, 2},
        true);
  }
  void initBox(HighFive::Group& _group) const;
  void initConnectivity(HighFive::File& _file,
                        const networkV4::network& _net) const;
  void initBondType(HighFive::File& _file,
                    const std::string& _name,
                    const networkV4::network& _net,
                    const networkV4::bondType _type) const;
  void initObservations(HighFive::File& _file,
                        const networkV4::network& _net) const;

private:
  void saveBox(HighFive::Group& _group,
               const networkV4::network& _net,
               const std::size_t _step,
               const double _time) const;
  void savePositions(HighFive::Group& _group,
                     const networkV4::network& _net,
                     const std::size_t _step,
                     const double _time) const;
  void saveConnected(HighFive::Group& _group,
                     const networkV4::network& _net,
                     const std::size_t _step,
                     const double _time,
                     const networkV4::bondType _type) const;

private:
  void saveConnectivity2Dataset(HighFive::DataSet& _dataset,
                                const networkV4::bonds& _bonds,
                                const networkV4::bondType _type) const;
  size_t findChunkSize(const size_t _size, const size_t _maxSize) const;

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
    _dataset.resize({_dataset.getSpace().getDimensions()[0] + 1});
    _dataset.select({_dataset.getSpace().getDimensions()[0] - 1}, {1})
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
                         const std::pair<size_t, size_t>& _size) const
      -> HighFive::DataSet
  {
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