#pragma once

#include <filesystem>

#include <highfive/highfive.hpp>

#include "Core/Network.hpp"
#include "NetworkOut.hpp"

class networkOutHDF5 : public networkOut
{
public:
  networkOutHDF5(const std::filesystem::path& _path,
                 const networkV4::network& _net);
  ~networkOutHDF5();

public:
  void save(const networkV4::network& _net,
            const std::size_t _step,
            const double _time,
            const std::string& _type) const override;

private:
  void writeHeader(HighFive::File& _file) const;
  void initParticles(HighFive::File& _file, const networkV4::network& _net);
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
  void saveScalar(HighFive::DataSet& _dataset, const T& _data) const
  {
    _dataset.resize({_dataset.getSpace().getDimensions()[0] + 1});
    _dataset.select({_dataset.getSpace().getDimensions()[0] - 1}, {1})
        .write(_data);
  }

  template<typename T>
  HighFive::DataSet createDataSet(HighFive::Group& _group,
                                  const std::string& _name,
                                  std::vector<std::size_t> _dims,
                                  std::vector<std::size_t> _maxDims,
                                  std::vector<hsize_t> _chunkSize,
                                  bool _deflate) const
  {
    HighFive::DataSetCreateProps props;
    props.add(HighFive::Chunking(_chunkSize));
    if (_deflate)
      props.add(HighFive::Deflate(9));
    return _group.createDataSet<T>(
        _name, HighFive::DataSpace(_dims, _maxDims), props);
  }
};