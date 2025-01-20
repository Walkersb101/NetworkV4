#pragma once

#include <cstddef>

#include "Integration/Integrators/Overdamped/Euler.hpp"
#include "Integration/Integrators/Overdamped/EulerHeun.hpp"
#include "Integration/Integrators/Overdamped/AdaptiveEulerHeun.hpp"

namespace networkV4
{
namespace integration
{

template<typename Integrator>
class Run
{
public:
  Run() = delete;
  Run(Integrator _integrator)
      : m_integrator(_integrator)
  {
  }
  ~Run() = default;

  template<typename StopFunction>
  void run(network& _network, std::size_t _steps, StopFunction _stopFunction)
  {
    size_t startStep = m_steps;
    size_t targetSteps = _steps + startStep;
    while (m_steps++ < targetSteps) {
      m_integrator.step(_network);
      m_t += m_integrator.getDt();
      if (_stopFunction(m_t, m_integrator.getDt(), _network))
        break;
    }
  }

private:
  Integrator m_integrator;
  double m_t = 0.0;
  size_t m_steps = 0;
};

}  // namespace Integration
}  // namespace networkV4
