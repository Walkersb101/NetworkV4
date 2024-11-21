#pragma once

#include "Integration/Integrator.hpp"
#include "Integration/OverDamped.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"
#include "Core/BreakTypes.hpp"

namespace networkV4
{

class SteepestDescent : public integrator
{
public:
  SteepestDescent();
  SteepestDescent(double _tol);
  ~SteepestDescent();

public:
  void integrate(network& _network) override;
  auto getDt() const -> double;

private:
  auto eGradient(network& _network, const double _step = 1e-6) const -> double;
  auto lineSearch(network& _network) const -> double;

private:
  double m_dt = config::integrators::default_dt;
  double m_tol = config::integrators::default_dt;

  double m_dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double m_dtMax = config::integrators::adaptiveIntegrator::dtMax;
  double m_rootTol = config::rootMethods::targetTol;

  size_t m_maxIter = config::integrators::miminizer::maxIter;
};

class FireMinimizer : public integrator
{
public:
  FireMinimizer();
  FireMinimizer(double _tol);
  FireMinimizer(double _dt, double _tol);
  ~FireMinimizer();

public:
  void integrate(network& _network) override;

private:
  auto xdoty(const std::vector<vec2d>& _x,
             const std::vector<vec2d>& _y) -> double;

private:
  double m_dt = config::integrators::default_dt;
  double m_tol = config::integrators::miminizer::tol;

private:
  size_t m_maxIter = config::integrators::miminizer::maxIter;
  double m_dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double m_dtMax = config::integrators::adaptiveIntegrator::dtMax;

  double m_alpha0 = config::integrators::fire2::alpha0;
  double m_finc = config::integrators::fire2::finc;
  double m_fdec = config::integrators::fire2::fdec;
  double m_falpha = config::integrators::fire2::falpha;
  size_t m_Ndelay = config::integrators::fire2::Ndelay;
  size_t m_Nnegmax = config::integrators::fire2::Nnegmax;
  double m_dmax = config::integrators::fire2::dmax;
};

class OverdampedAdaptiveMinimizer : public integrator
{
public:
  OverdampedAdaptiveMinimizer();
  OverdampedAdaptiveMinimizer(double _esp, double _tol);
  ~OverdampedAdaptiveMinimizer();

public:
  void integrate(network& _network) override;
  auto getDt() const -> double;

public:
  auto innerIteration(network& _network, double _dt) -> double;

private:
  double m_dt = config::integrators::default_dt;
  double m_esp = config::integrators::adaptiveIntegrator::esp;
  double m_tol = config::integrators::miminizer::tol;

  size_t m_maxInnerIter = config::integrators::adaptiveIntegrator::maxIter;
  size_t m_maxOuterIter = config::integrators::miminizer::maxIter;

  std::vector<vec2d> m_frn;
  std::vector<vec2d> m_frnbar;

  double m_qMin = config::integrators::adaptiveHeun::qMin;
  double m_qMax = config::integrators::adaptiveHeun::qMax;
  double m_dtMin = config::integrators::adaptiveIntegrator::dtMin;
  double m_dtMax = config::integrators::adaptiveIntegrator::dtMax;
  double m_fdec = config::integrators::OverdampedAdaptiveMinimizer::fdec;
  double m_dmax = config::integrators::OverdampedAdaptiveMinimizer::dmax;
};

}  // namespace networkV4