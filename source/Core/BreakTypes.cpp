#include <iostream>
#include <random>

#include "BreakTypes.hpp"

networkV4::BreakTypes::BreakTypes() {}

networkV4::BreakTypes::~BreakTypes() {}

networkV4::NoBreaks::NoBreaks() {}

networkV4::NoBreaks::~NoBreaks() {}

void networkV4::NoBreaks::Data(const network& _network,
                               const integrator& _integrator,
                               double& _maxVal,
                               std::size_t& _count)
{
  _maxVal = 0.0;
  _count = 0;
}

auto networkV4::NoBreaks::Break(
    network& _network, const integrator& _integrator) -> std::vector<bond*>
{
  return std::vector<bond*>();
}

networkV4::StrainBreak::StrainBreak() {}

networkV4::StrainBreak::~StrainBreak() {}

void networkV4::StrainBreak::Data(const network& _network,
                                  const integrator& _integrator,
                                  double& _maxVal,
                                  std::size_t& _count)
{
  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (const auto& b : connectedList) {
    const double strain = _network.bondStrain(b);
    _maxVal = std::max(_maxVal, strain - b.lambda());
    if (strain >= b.lambda()) {
      ++_count;
    }
  }
}

auto networkV4::StrainBreak::Break(
    network& _network, const integrator& _integrator) -> std::vector<bond*>
{
  std::vector<bond*> broken;
  auto& bonds = _network.getBonds();
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (auto& b : connectedList) {
    const double strain = _network.bondStrain(b);
    if (strain < b.lambda()) {
      b.connected() = false;
      broken.push_back(&b);
    }
  }
  return broken;
}

networkV4::EnergyBreak::EnergyBreak() {}

networkV4::EnergyBreak::~EnergyBreak() {}

void networkV4::EnergyBreak::Data(const network& _network,
                                  const integrator& _integrator,
                                  double& _maxVal,
                                  std::size_t& _count)
{
  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (const auto& b : connectedList) {
    const double strain = _network.bondEnergy(b);
    _maxVal = std::max(_maxVal, strain - b.lambda());
    if (strain >= b.lambda()) {
      ++_count;
    }
  }
}

auto networkV4::EnergyBreak::Break(
    network& _network, const integrator& _integrator) -> std::vector<bond*>
{
  std::vector<bond*> broken;
  auto& bonds = _network.getBonds();
  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (auto& b : connectedList) {
    const double strain = _network.bondEnergy(b);
    if (strain < b.lambda()) {
      b.connected() = false;
      broken.push_back(&b);
    }
  }
  return broken;
}

networkV4::SGRBreak::SGRBreak() {}

networkV4::SGRBreak::SGRBreak(double _T, double _tau0, unsigned long _seed)
    : m_T(_T)
    , m_rT(1.0 / _T)
    , m_tau0(_tau0)
    , m_rtau0(1.0 / _tau0)
    , m_rng(_seed)
{
}

networkV4::SGRBreak::~SGRBreak() {}

void networkV4::SGRBreak::Data(const network& _network,
                               const integrator& _integrator,
                               double& _maxVal,
                               std::size_t& _count)
{
  auto rng = m_rng;

  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;

  std::uniform_real_distribution<double> dist(0.0, 1.0);

  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (const auto& b : connectedList) {
    const double energy = _network.bondEnergy(b);
    const double distAbove = energy - b.lambda();
    _maxVal = std::max(_maxVal, distAbove);

    const double rate = m_rtau0 * std::exp(-distAbove * m_rT);
    const double probOfBreak = 1 - std::exp(-rate * _integrator.getDt());
    const double rand = dist(m_rng);

    if (rand < probOfBreak) {
      _count++;
    }
  }
}

auto networkV4::SGRBreak::Break(network& _network,
                                const integrator& _integrator)
    -> std::vector<bond*>
{
  auto& bonds = _network.getBonds();
  std::vector<bond*> broken;
  std::uniform_real_distribution<double> dist(0.0, 1.0);

  auto connected = [](const auto& _b) { return _b.connected(); };
  auto connectedList = bonds | std::views::filter(connected);
  for (auto& b : connectedList) {
    const double energy = _network.bondEnergy(b);
    const double distAbove = energy - b.lambda();

    const double rate = m_rtau0 * std::exp(-distAbove * m_rT);
    const double probOfBreak = 1 - std::exp(-rate * _integrator.getDt());
    const double rand = dist(m_rng);

    if (rand < probOfBreak) {
      b.connected() = false;
      broken.push_back(&b);
    }
  }
  return broken;
}