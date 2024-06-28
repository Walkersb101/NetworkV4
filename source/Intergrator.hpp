#pragma once

#include "Network.hpp"
#include "Nodes.hpp"

namespace network
{
class intergrator
{
public:
  intergrator();

  virtual ~intergrator();

public:
  virtual void intergrate(network& _network);
};

class overdampedEuler : public intergrator
{
public:
  overdampedEuler();
  overdampedEuler(double _dt);
  ~overdampedEuler();

public:
  void intergrate(network& _network) override;

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
  void intergrate(network& _network) override;

private:
  double m_dt;
  std::vector<vec2d> m_tempForces;
  std::unordered_map<bondType, tensor2d> m_tempStresses;
};

}  // namespace network
