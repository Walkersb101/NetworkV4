#pragma once

#include <optional>

#include "Misc/Vec2.hpp"

namespace networkV4
{
namespace BreakTypes
{
class StrainBreak
{
public:
  StrainBreak(double _lambda, double _r0)
      : m_lambda(_lambda)
      , m_r0(_r0)
  {
    if (_r0 <= 0.0) {
      throw("StrainBreak: r0 must be positive");
    }
    m_invR0 = 1.0 / m_r0;
  }

public:
  const double lambda() const { return m_lambda; }
  const double r0() const { return m_r0; }

public:
  bool checkBreak(const Utils::vec2d& _r) const
  {
    return (_r.norm() * m_invR0) - 1.0 > m_lambda;
  }
  std::optional<double> thresholdData(const Utils::vec2d& _r) const
  {
    return (_r.norm() * m_invR0) - 1.0 - m_lambda;
  }
  std::optional<double> data(const Utils::vec2d& _r) const
  {
    return (_r.norm() * m_invR0) - 1.0;
  }

private:
  double m_r0;  // equilibrium bond length
  double m_invR0;
  double m_lambda;
};
}  // namespace BreakTypes
}  // namespace networkV4