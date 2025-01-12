#pragma once

#include <optional>

#include "Misc/Config.hpp"
#include "Misc/Math/Vector.hpp"

namespace networkV4
{
namespace Forces
{
class HarmonicBond
{
public:
  HarmonicBond(double _k, double _r0, bool _normalized = false)
      : m_k(_k)
      , m_r0(_r0)
      , m_normalized(_normalized)
  {
    if (m_normalized) {
      if (m_r0 < ROUND_ERROR_PRECISION) {
        throw("HarmonicBond: r0 is too small");
      }
      m_k /= m_r0;
    }
  }

public:
  const double k() const { return m_normalized ? m_k * m_r0 : m_k; }
  const double r0() const { return m_r0; }
  const bool normalized() const { return m_normalized; }

public:
  std::optional<Utils::Math::vec2d> force(const Utils::Math::vec2d& _dx) const
  {
    const auto r = _dx.norm();
    const auto dr = r - m_r0;
    auto fac = -m_k * dr;
    if (r > ROUND_ERROR_PRECISION) {
      fac /= r;
    } else {
      throw("HarmonicBond::force: r is too small");
    }
    return _dx * fac;
  }
  std::optional<double> energy(const Utils::Math::vec2d& _dx) const
  {
    const auto r = _dx.norm();
    const auto dr = r - m_r0;
    return 0.5 * m_k * dr * dr;
  }

private:
  double m_k;  // spring constant
  double m_r0;  // equilibrium bond length
  bool m_normalized = false;  // whether the spring constant is normalized
};
}  // namespace Forces
}  // namespace networkV4