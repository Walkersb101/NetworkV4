#pragma once

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

namespace networkV4
{
namespace integration
{

class IntegratorBase
{
public:
  IntegratorBase() = default;
  IntegratorBase(double _dt)
      : m_dt(_dt)
  {
  }
  virtual ~IntegratorBase() = default;

public:
  virtual void step(network& _network) = 0;
  auto getDt() const -> double { return m_dt; }

public:
  double m_dt = config::integrators::default_dt;
};

}  // namespace integration
}  // namespace networkV4
