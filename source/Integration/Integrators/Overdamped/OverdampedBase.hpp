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
  Overdamped(double _zeta)
      : m_zeta(_zeta)
      , m_invZeta(1.0 / _zeta)
  {
  }
  virtual ~Overdamped() = default;

public:
    auto getZeta() const -> double { return m_zeta; }

protected:
  double m_zeta = 1.0;
  double m_invZeta = 1.0;
};

}  // namespace integration
}  // namespace networkV4