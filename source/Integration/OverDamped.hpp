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
  double m_dt = config::integrators::default_dt;
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
  double m_dt = config::integrators::default_dt;
  std::vector<vec2d> m_frn;
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

private:
  double m_dt = config::integrators::default_dt;
  double m_nextdt = config::integrators::default_dt;

  double m_esp = config::integrators::adaptiveIntegrator::esp;

  std::vector<vec2d> m_frn;
  std::vector<vec2d> m_frnbar;

  double m_qMin = config::integrators::adaptiveHeun::qMin;
  double m_qMax = config::integrators::adaptiveHeun::qMax;
  double m_dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double m_dtMax = config::integrators::adaptiveIntegrator::dtMax;
  size_t m_maxInnerIter = config::integrators::adaptiveIntegrator::maxIter;
};

}  // namespace networkV4