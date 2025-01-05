#include <fstream>
#include <iostream>
#include <ranges>
#include <stdexcept>
#include <vector>

#include "hdf5Out.hpp"

#include <highfive/highfive.hpp>

void networkOutHDF5::initParticles(HighFive::File _file,
                                   const networkV4::network& _net)
{
  
}

void networkOutHDF5::initBox(HighFive::Group& _group) const
{
  auto box = createOrGetGroup(_group, "box");
  box.createAttribute("dimension", 2);
  box.createAttribute("boundary", std::array<bool, 2> {true, true});

  auto edges = createOrGetGroup(box, "edges");
  initScalarDataset<int>(edges, "step");
  initScalarDataset<double>(edges, "time");
  initMatrixDataset<double>(edges, "value");
}

void networkOutHDF5::initConnectivity(HighFive::File& _file,
                                      const networkV4::network& _net) const
{
  initBondType(_file, "sacrificial", _net, networkV4::bondType::sacrificial);
  initBondType(_file, "matrix", _net, networkV4::bondType::matrix);
}

void networkOutHDF5::initBondType(HighFive::File& _file,
                                  const std::string& _name,
                                  const networkV4::network& _net,
                                  const networkV4::bondType _type) const
{
  const size_t count = std::count_if(_net.getBonds().begin(),
                                     _net.getBonds().end(),
                                     [&_type](const networkV4::bond& b)
                                     { return b.type() == _type; });

  const size_t posChunkSize = std::min(config::hdf5::maxChunk, count);

  auto isvalid = [&_type](const auto& _b) { return _b.type() == _type; };
  auto connectedList = _net.getBonds() | std::views::filter(isvalid);

  // Save connectivity

  HighFive::Group connectivity =
      createOrGetGroup(_file, "particles/all/connectivity");

  auto value = createDataSet<int>(connectivity,
                                  _name,
                                  std::vector<std::size_t> {count, 2},
                                  std::vector<std::size_t> {count, 2},
                                  std::vector<hsize_t> {count, 2},
                                  true);
  value.createAttribute(
      "particles_group",
      HighFive::Reference(createOrGetGroup(_file, "particles"),
                          createOrGetGroup(_file, "particles/all")));

  std::vector<std::vector<size_t>> connections;
  connections.reserve(count);
  for (const networkV4::bond& bond : connectedList) {
    connections.push_back({bond.src(), bond.dst()});
  }
  value.write(connections);

  // Save bond parameters

  auto bondParameters = createOrGetGroup(_file, "parameters/bonds");
  auto typeGroup = createOrGetGroup(bondParameters, _name);

  auto mu = createDataSet<double>(typeGroup,
                                  "mu",
                                  std::vector<std::size_t> {count},
                                  std::vector<std::size_t> {count},
                                  std::vector<hsize_t> {count},
                                  true);
  auto r0 = createDataSet<double>(typeGroup,
                                  "r0",
                                  std::vector<std::size_t> {count},
                                  std::vector<std::size_t> {count},
                                  std::vector<hsize_t> {count},
                                  true);
  auto lambda = createDataSet<double>(typeGroup,
                                      "lambda",
                                      std::vector<std::size_t> {count},
                                      std::vector<std::size_t> {count},
                                      std::vector<hsize_t> {count},
                                      true);
  std::vector<double> data(count, 0.0);
  std::transform(connectedList.begin(),
                 connectedList.end(),
                 data.begin(),
                 [](const networkV4::bond& _bond) { return _bond.mu(); });
  mu.write(data);
  std::transform(connectedList.begin(),
                 connectedList.end(),
                 data.begin(),
                 [](const networkV4::bond& _bond)
                 { return _bond.naturalLength(); });
  r0.write(data);
  std::transform(connectedList.begin(),
                 connectedList.end(),
                 data.begin(),
                 [](const networkV4::bond& _bond) { return _bond.lambda(); });
  lambda.write(data);
}

void networkOutHDF5::initObservations(HighFive::File& _file,
                                      const networkV4::network& _net) const
{
  const size_t matrixCount =
      std::count_if(_net.getBonds().begin(),
                    _net.getBonds().end(),
                    [](const networkV4::bond& b)
                    { return b.type() == networkV4::bondType::matrix; });
  const size_t sacrificialCount =
      std::count_if(_net.getBonds().begin(),
                    _net.getBonds().end(),
                    [](const networkV4::bond& b)
                    { return b.type() == networkV4::bondType::sacrificial; });

  const size_t matrixChunkSize =
      findChunkSize(matrixCount, config::hdf5::maxChunk);
  const size_t sacrificialChunkSize =
      findChunkSize(sacrificialCount, config::hdf5::maxChunk);

  // Connected Observables
  auto sacrificialConnectedObs =
      createOrGetGroup(_file, "observables/sacrificialConnected");
  createDataSet<int>(sacrificialConnectedObs,
                     "step",
                     std::vector<std::size_t> {0},
                     std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
                     std::vector<hsize_t> {1},
                     false);
  createDataSet<double>(
      sacrificialConnectedObs,
      "time",
      std::vector<std::size_t> {0},
      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
      std::vector<hsize_t> {1},
      false);
  createDataSet<bool>(sacrificialConnectedObs,
                      "value",
                      std::vector<std::size_t> {0, sacrificialCount},
                      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED,
                                                sacrificialCount},
                      std::vector<hsize_t> {1, sacrificialChunkSize},
                      true);

  auto matrixConnectedObs =
      createOrGetGroup(_file, "observables/matrixConnected");
  createDataSet<int>(matrixConnectedObs,
                     "step",
                     std::vector<std::size_t> {0},
                     std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
                     std::vector<hsize_t> {1},
                     false);
  createDataSet<double>(
      matrixConnectedObs,
      "time",
      std::vector<std::size_t> {0},
      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
      std::vector<hsize_t> {1},
      false);
  createDataSet<bool>(
      matrixConnectedObs,
      "value",
      std::vector<std::size_t> {0, matrixCount},
      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED, matrixCount},
      std::vector<hsize_t> {1, matrixChunkSize},
      true);

  // type Observables
  auto typeObs = createOrGetGroup(_file, "observables/logType");
  createDataSet<int>(typeObs,
                     "step",
                     std::vector<std::size_t> {0},
                     std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
                     std::vector<hsize_t> {1},
                     false);
  createDataSet<double>(
      typeObs,
      "time",
      std::vector<std::size_t> {0},
      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
      std::vector<hsize_t> {1},
      false);
  createDataSet<std::string>(
      typeObs,
      "value",
      std::vector<std::size_t> {0},
      std::vector<std::size_t> {HighFive::DataSpace::UNLIMITED},
      std::vector<hsize_t> {1},
      false);
}

void networkOutHDF5::saveBox(HighFive::Group& _group,
                             const networkV4::network& _net,
                             const std::size_t _step,
                             const double _time) const
{
  auto step = _group.getDataSet("step");
  saveScalar(step, _step);

  auto time = _group.getDataSet("time");
  saveScalar(time, _time);

  auto value = _group.getDataSet("value");
  const std::vector<std::vector<double>> data = {
      {_net.getDomain().x, 0},
      {_net.getDomain().y * _net.getShearStrain(), _net.getDomain().y}};
  value.resize({time.getSpace().getDimensions()[0], 2, 2});
  value.select({time.getSpace().getDimensions()[0] - 1, 0, 0}, {1, 2, 2})
      .squeezeMemSpace({0})
      .write(data);
}

void networkOutHDF5::savePositions(HighFive::Group _group,
                                   const networkV4::network& _net,
                                   const std::size_t _step,
                                   const double _time) const
{
  auto step = _group.getDataSet("step");
  saveScalar(step, _step);

  auto time = _group.getDataSet("time");
  saveScalar(time, _time);

  const double offset = _net.getShearStrain() * _net.getDomain().y * 0.5;
  std::vector<std::vector<double>> data;
  data.reserve(_net.getNodes().size());
  for (const auto& pos : _net.getNodes().positions()) {
    const vec2d shearOffset =
        vec2d((_net.getShearStrain() * pos.y) - offset, 0.0);
    const vec2d TriPos = _net.wrapedPosition(pos - shearOffset) + shearOffset;
    data.push_back({TriPos.x, TriPos.y});
  }

  auto value = _group.getDataSet("value");
  // const double* data =
  //     reinterpret_cast<const double*>(_net.getNodes().positions().data());

  value.resize({time.getSpace().getDimensions()[0], _net.getNodes().size(), 2});
  // value
  //     .select({time.getSpace().getDimensions()[0] - 1, 0, 0},
  //             {1, _net.getNodes().size(), 2})
  //     .squeezeMemSpace({0})
  //     .write_raw(data);
  value
      .select({time.getSpace().getDimensions()[0] - 1, 0, 0},
              {1, _net.getNodes().size(), 2})
      .squeezeMemSpace({0})
      .write(data);
}

void networkOutHDF5::saveConnected(HighFive::Group& _group,
                                   const networkV4::network& _net,
                                   const std::size_t _step,
                                   const double _time,
                                   const networkV4::bondType _type) const
{
  const size_t count = std::count_if(_net.getBonds().begin(),
                                     _net.getBonds().end(),
                                     [&_type](const networkV4::bond& b)
                                     { return b.type() == _type; });

  auto isvalid = [&_type](const auto& _b) { return _b.type() == _type; };
  auto connectedList = _net.getBonds() | std::views::filter(isvalid);

  auto step = _group.getDataSet("step");
  saveScalar(step, _step);

  auto time = _group.getDataSet("time");
  saveScalar(time, _time);

  std::vector<bool> data(count, false);
  std::transform(connectedList.begin(),
                 connectedList.end(),
                 data.begin(),
                 [](const networkV4::bond& _bond)
                 { return _bond.connected(); });

  auto value = _group.getDataSet("value");
  value.resize({time.getSpace().getDimensions()[0], count});
  value.select({time.getSpace().getDimensions()[0] - 1, 0}, {1, count})
      .squeezeMemSpace({0})
      .write(data);
}

size_t networkOutHDF5::findChunkSize(const size_t _size,
                                     const size_t _maxSize) const
{
  size_t maxi = 1;
  if (_size < _maxSize)
    return _size;
  for (size_t i = 2; i <= _maxSize; i++) {
    if (_size % i == 0)
      maxi = std::max(maxi, i);
  }
  return maxi;
}