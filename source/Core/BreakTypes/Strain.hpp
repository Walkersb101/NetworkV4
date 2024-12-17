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
  StrainBreak(double _lambda)     : m_lambda(_lambda)
{
}

public:
  const double lambda() const { return m_lambda; }

public:
  bool checkBreak(const Utils::vec2d& _r) const {
  return _r.norm() > m_lambda;
}
  std::optional<double> thresholdData(const Utils::vec2d& _r) const {
  return _r.norm() - m_lambda;
}

private:
  double m_lambda;
};
}  // namespace BreakTypes
}  // namespace networkV4