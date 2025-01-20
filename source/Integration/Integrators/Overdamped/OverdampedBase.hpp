#pragma once

#pragma once

#include <algorithm>
#include <vector>

namespace networkV4
{
namespace integration
{

class Overdamped
{
public:
  Overdamped() = default;
  Overdamped(double _gamma)
      : m_gamma(_gamma)
      , m_invGamma(1.0 / _gamma)
  {
  }
  virtual ~Overdamped() = default;

protected:
  double m_gamma = 1.0;
  double m_invGamma = 1.0;
};

}  // namespace integration
}  // namespace networkV4