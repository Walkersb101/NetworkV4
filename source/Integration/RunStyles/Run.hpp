#pragma once

#include <cstddef>

namespace networkV4
{
namespace Integration
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
    for (std::size_t i = 0; i < _steps; ++i) {
      m_integrator.step(_network);
      m_t += m_integrator.getDt();
      if (_stopFunction(m_t, m_integrator.getDt(), _network))
        break;
    }
  }

private:
  Integrator m_integrator;
  double m_t = 0.0;
};

}  // namespace Integration
}  // namespace networkV4
