#pragma once

#include "Network.hpp"
#include "Nodes.hpp"
#include "Tensor2.hpp"

namespace networkV4
{
    enum class intergratorType : std::uint8_t
    {
        overdampedEuler,
        overdampedEulerHeun,
        overdampedAdaptiveEulerHeun
    };

class intergrator
{
public:
  intergrator();

  virtual ~intergrator();

public:
  virtual void integrate(network& _network);
};

class overdampedEuler : public intergrator
{
public:
  overdampedEuler();
  overdampedEuler(double _dt);
  ~overdampedEuler();

public:
  void integrate(network& _network) override;

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

}  // namespace networkV4
