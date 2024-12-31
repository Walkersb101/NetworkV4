#pragma once

#include <cstddef>

namespace networkV4
{
namespace Integration
{

template<typename integrator>
class Run
{
public:
  Run() = default;
  ~Run() = default;
  void run(std::size_t _steps, auto _stopFunction)
  {
    for (std::size_t i = 0; i < _steps; ++i) {
      integrator.step();
      if (_stopFunction())
        break;
    }
  }

private: 
integrator m_integrator;
double m_t;

};

}  // namespace Integration
}  // namespace networkV4
