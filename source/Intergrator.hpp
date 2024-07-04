#pragma once

#include <cstddef>
#include <cstdint>
#include <unordered_map>
#include <vector>

#include "Network.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"
#include "Vec2.hpp"

namespace networkV4
{
enum class intergratorType : std::uint8_t
{
  overdampedEuler,
  overdampedEulerHeun,
  overdampedAdaptiveEulerHeun,
  FireMinimizer,
};

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
  std::unordered_map<bondType, tensor2d> m_tempStresses;
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
  std::unordered_map<bondType, tensor2d> m_tempStresses;
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
  auto power(const network& _network) -> double;

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
