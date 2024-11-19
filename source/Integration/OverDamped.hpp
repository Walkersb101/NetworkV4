#pragma once

#include "Integration/Integrator.hpp"

namespace networkV4
{

class overdampedEuler : public integrator
{
public:
  overdampedEuler();
  overdampedEuler(double _dt);
  ~overdampedEuler();

public:
  void integrate(network& _network) override;
  auto getDt() const -> double;

private:
  double m_dt;
};

class OverdampedEulerHeun : public integrator
{
public:
  OverdampedEulerHeun();
  OverdampedEulerHeun(double _dt);
  ~OverdampedEulerHeun();

public:
  void integrate(network& _network) override;
  auto getDt() const -> double;

private:
  double m_dt;
  std::vector<vec2d> m_tempPositions;
  std::vector<vec2d> m_tempForces;
  stressMap m_tempStresses;
};

class OverdampedAdaptiveEulerHeun : public integrator
{
public:
  OverdampedAdaptiveEulerHeun();
  OverdampedAdaptiveEulerHeun(double _dt, double _esp);
  ~OverdampedAdaptiveEulerHeun();

public:
  void integrate(network& _network) override;
  auto getDt() const -> double;

private:
  auto forceErrorNorm(nodes& _network) -> double;

private:
  double m_dt;
  double m_nextdt;

  double m_esp;

  std::vector<vec2d> m_tempForces;
  std::vector<vec2d> m_tempPositions;
  stressMap m_tempStresses;
};

}  // namespace networkV4