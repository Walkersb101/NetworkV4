#pragma once

#include "Integration/Integrator.hpp"
#include "Integration/OverDamped.hpp"
#include "Misc/Roots.hpp"
#include "Misc/Tools.hpp"

namespace networkV4
{

class SteepestDescent : public integrator
{
public:
  SteepestDescent();
  SteepestDescent(double _dt);
  ~SteepestDescent();

public:
    void integrate(network& _network) override;
    auto getDt() const -> double;

private:
    auto eGradient(network& _network, const double _step = 1e-6) const -> double;
    auto lineSearch(network& _network) const -> double;

private:
    double m_dt;
    double m_nextdt;
    double m_prevEnergy;
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
  double m_dt;
  double m_tol;
};

class OverdampedAdaptiveMinimizer : public integrator
{
public:
  OverdampedAdaptiveMinimizer();
  OverdampedAdaptiveMinimizer(double _dt, double _esp, double _tol);
  ~OverdampedAdaptiveMinimizer();

public:
  void integrate(network& _network) override;

private:
  OverdampedAdaptiveEulerHeun m_integrator;

  double m_tol;
};

}  // namespace networkV4