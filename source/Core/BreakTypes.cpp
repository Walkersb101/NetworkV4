#include "BreakTypes.hpp"

#include <random>
#include <iostream>

networkV4::BreakTypes::BreakTypes() {}

networkV4::BreakTypes::~BreakTypes() {}

networkV4::NoBreaks::NoBreaks() {}

networkV4::NoBreaks::~NoBreaks() {}

void networkV4::NoBreaks::Data(const network& _network,
                               const intergrator& _intergrator,
                               double& _maxVal,
                               std::size_t& _count)
{
  _maxVal = 0.0;
  _count = 0;
}

auto networkV4::NoBreaks::Break(network& _network,
                                const intergrator& _intergrator)
    -> std::vector<std::size_t>
{
  return std::vector<std::size_t>();
}

networkV4::StrainBreak::StrainBreak()
{
  throw std::runtime_error("lambda needs to be specified ");
}

networkV4::StrainBreak::StrainBreak(const std::vector<double>& _lambda)
    : m_lambda(_lambda)
{
}

networkV4::StrainBreak::~StrainBreak() {}

void networkV4::StrainBreak::Data(const network& _network,
                                  const intergrator& _intergrator,
                                  double& _maxVal,
                                  std::size_t& _count)
{
  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;
  for (size_t i = 0; i < bonds.size(); ++i) {
    const auto& b = bonds.get(i);
    if (!b.connected()) {
      continue;
    }
    const double strain = _network.bondStrain(b);
    const double distAbove = strain - m_lambda[i];
    _maxVal = std::max(_maxVal, distAbove);
    if (distAbove >= 0.0) {
      ++_count;
    }
  }
}

auto networkV4::StrainBreak::Break(network& _network,
                                   const intergrator& _intergrator)
    -> std::vector<std::size_t>
{
  auto& bonds = _network.getBonds();
  std::vector<std::size_t> broken;
  for (std::size_t i = 0; i < bonds.size(); ++i) {
    auto& b = bonds.get(i);
    if (!b.connected() || _network.bondStrain(b) < m_lambda[i]) {
      continue;
    }
    b.connected() = false;
    broken.push_back(i);
  }
  return broken;
}

networkV4::EnergyBreak::EnergyBreak()
{
  throw std::runtime_error("Energy needs to be specified ");
}

networkV4::EnergyBreak::EnergyBreak(const std::vector<double>& _E)
    : m_E(_E)
{
}

networkV4::EnergyBreak::~EnergyBreak() {}

void networkV4::EnergyBreak::Data(const network& _network,
                                  const intergrator& _intergrator,
                                  double& _maxVal,
                                  std::size_t& _count)
{
  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;
  for (size_t i = 0; i < bonds.size(); ++i) {
    const auto& b = bonds.get(i);
    if (!b.connected()) {
      continue;
    }
    const double energy = _network.bondEnergy(b);
    const double distAbove = energy - m_E[i];
    _maxVal = std::max(_maxVal, distAbove);
    if (distAbove >= 0.0) {
      ++_count;
    }
  }
}

auto networkV4::EnergyBreak::Break(network& _network,
                                   const intergrator& _intergrator)
    -> std::vector<std::size_t>
{
  auto& bonds = _network.getBonds();
  std::vector<std::size_t> broken;
  for (std::size_t i = 0; i < bonds.size(); ++i) {
    auto& b = bonds.get(i);
    if (!b.connected() || _network.bondEnergy(b) < m_E[i]) {
      continue;
    }
    b.connected() = false;
    broken.push_back(i);
  }
  return broken;
}

networkV4::SGRBreak::SGRBreak() {}

networkV4::SGRBreak::SGRBreak(const std::vector<double>& _E,
                              double _T,
                              double _tau0,
                              unsigned long _seed)
    : m_E(_E)
    , m_T(_T)
    , m_rT(1.0 / _T)
    , m_tau0(_tau0)
    , m_rtau0(1.0 / _tau0)
    , m_rng(_seed)
{
}

networkV4::SGRBreak::~SGRBreak() {}

void networkV4::SGRBreak::Data(const network& _network,
                               const intergrator& _intergrator,
                               double& _maxVal,
                               std::size_t& _count)
{
  const auto& bonds = _network.getBonds();
  _count = 0;
  _maxVal = -1e10;
  auto rng = m_rng;
  for (size_t i = 0; i < bonds.size(); i++) {
    const auto& b = bonds.get(i);
    if (!b.connected()) {
      continue;
    }
    const double energy = _network.bondEnergy(b);
    const double distAbove = energy - m_E[i];
    const double rate = m_rtau0 * std::exp(-distAbove * m_rT);
    std::exponential_distribution<double> dist(rate);
    if (dist(rng) >= _intergrator.getDt()) {
      ++_count;
    }
    _maxVal = std::max(_maxVal, distAbove);
  }
}

auto networkV4::SGRBreak::Break(network& _network,
                                const intergrator& _intergrator)
    -> std::vector<std::size_t>
{
  auto& bonds = _network.getBonds();
  std::vector<std::size_t> broken;
  std::uniform_real_distribution<double> dist(0.0, 1.0);
  for (size_t i = 0; i < bonds.size(); i++) {
    auto& b = bonds.get(i);
    if (!b.connected()) {
      continue;
    }
    const double energy = _network.bondEnergy(b);
    const double distAbove = energy - m_E[i];
    const double rate = m_rtau0 * std::exp(-distAbove * m_rT);
    const double probOfBreak = 1 - std::exp(-rate * _intergrator.getDt());
    const double rand = dist(m_rng);
    if (rand < probOfBreak) {
      b.connected() = false;
      broken.push_back(i);
    }
  }
    return broken;
}