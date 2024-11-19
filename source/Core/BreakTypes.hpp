#pragma once

#include <trng/yarn2.hpp>

#include "Misc/Enums.hpp"
#include "Intergration/Intergrator.hpp"
#include "Core/Network.hpp"

namespace networkV4
{
class BreakTypes
{
public:
  BreakTypes();
  virtual ~BreakTypes();

public:
  virtual void Data(const network& _network,
                    const intergrator& _intergrator,
                    double& _maxVal,
                    std::size_t& _count) = 0;
  virtual auto Break(network& _network, const intergrator& _intergrator)
      -> std::vector<std::size_t> = 0;
};

class NoBreaks : public BreakTypes
{
public:
  NoBreaks();
  virtual ~NoBreaks();

public:
  void Data(const network& _network,
            const intergrator& _intergrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const intergrator& _intergrator)
      -> std::vector<std::size_t> override;
};

class StrainBreak : public BreakTypes
{
public:
  StrainBreak();
  StrainBreak(const std::vector<double>& _lambda);
  virtual ~StrainBreak();

public:
  void Data(const network& _network,
            const intergrator& _intergrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const intergrator& _intergrator)
      -> std::vector<std::size_t> override;

private:
  std::vector<double> m_lambda;
};

class EnergyBreak : public BreakTypes
{
public:
  EnergyBreak();
  EnergyBreak(const std::vector<double>& _E);
  virtual ~EnergyBreak();

public:
  void Data(const network& _network,
            const intergrator& _intergrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const intergrator& _intergrator)
      -> std::vector<std::size_t> override;

private:
  std::vector<double> m_E;
};

class SGRBreak : public BreakTypes
{
public:
  SGRBreak();
  SGRBreak(const std::vector<double>& _E, double _T, double _tau0, unsigned long _seed);
  virtual ~SGRBreak();

public:
  void Data(const network& _network,
            const intergrator& _intergrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const intergrator& _intergrator)
      -> std::vector<std::size_t> override;

private:
std::vector<double> m_E;
  double m_T;
  double m_rT;
  double m_tau0;
  double m_rtau0;

  trng::yarn2 m_rng;
};

}  // namespace networkV4