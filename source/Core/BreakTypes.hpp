#pragma once

#include <trng/yarn2.hpp>

#include "Core/Network.hpp"
#include "Integration/Integrator.hpp"
#include "Misc/Enums.hpp"

namespace networkV4
{
class BreakTypes
{
public:
  BreakTypes();
  virtual ~BreakTypes();

public:
  virtual void Data(const network& _network,
                    const integrator& _integrator,
                    double& _maxVal,
                    std::size_t& _count) = 0;
  virtual auto Break(network& _network, const integrator& _integrator)
      -> std::vector<bond*> = 0;
};

class NoBreaks : public BreakTypes
{
public:
  NoBreaks();
  virtual ~NoBreaks();

public:
  void Data(const network& _network,
            const integrator& _integrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const integrator& _integrator)
      -> std::vector<bond*> override;
};

class StrainBreak : public BreakTypes
{
public:
  StrainBreak();
  virtual ~StrainBreak();

public:
  void Data(const network& _network,
            const integrator& _integrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const integrator& _integrator)
      -> std::vector<bond*> override;
};

class EnergyBreak : public BreakTypes
{
public:
  EnergyBreak();
  virtual ~EnergyBreak();

public:
  void Data(const network& _network,
            const integrator& _integrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const integrator& _integrator)
      -> std::vector<bond*> override;
};

class SGRBreak : public BreakTypes
{
public:
  SGRBreak();
  SGRBreak(double _T,
           double _tau0,
           unsigned long _seed);
  virtual ~SGRBreak();

public:
  void Data(const network& _network,
            const integrator& _integrator,
            double& _maxVal,
            std::size_t& _count) override;
  auto Break(network& _network, const integrator& _integrator)
      -> std::vector<bond*> override;

private:
  double m_T;
  double m_rT;
  double m_tau0;
  double m_rtau0;

  trng::yarn2 m_rng;
};

}  // namespace networkV4