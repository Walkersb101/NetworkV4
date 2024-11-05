#pragma once

#include <cstddef>
#include <cstdint>
#include <vector>

#include "Network.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"
#include "EnumMap.hpp"

namespace networkV4
{

using stressMap = EnumMap<bondType, tensor2d, bondType::single, bondType::matrix>;

class intergrator
{
public:
  intergrator();

  virtual ~intergrator();

public:
  virtual void integrate(network& _network);
  virtual auto getDt() const -> double;
};

class overdampedEuler : public intergrator
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

class OverdampedEulerHeun : public intergrator
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

class OverdampedAdaptiveEulerHeun : public intergrator
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

class SteepestDescent : public intergrator
{
public:
  SteepestDescent();
  SteepestDescent(double _dt);
  ~SteepestDescent();

public:
    void integrate(network& _network) override;
    auto getDt() const -> double;

private:
    //auto eGradient(network& _network, const double _step = 0.01) const -> double;
    //auto lineSearch(network& _network) const -> double;

private:
    double m_dt;
    double m_nextdt;
    double m_prevEnergy;
};

class FireMinimizer : public intergrator
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

class OverdampedAdaptiveMinimizer : public intergrator
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
