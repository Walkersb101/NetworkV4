#pragma once

#include <algorithm>
#include <vector>

#include <range/v3/view/zip.hpp>

#include "Core/Network.hpp"
#include "Integration/Integrators/Adaptive.hpp"
#include "Integration/Integrators/IntegratorBase.hpp"
#include "Integration/Integrators/Overdamped/OverdampedBase.hpp"

namespace networkV4
{
namespace integration
{
// based on https://arxiv.org/pdf/2108.03399
class AdaptiveOverdampedEulerHeun
    : public Overdamped
    , public Adaptive
    , public IntegratorBase
{
public:
  AdaptiveOverdampedEulerHeun() = default;
  AdaptiveOverdampedEulerHeun(double _gamma,
                              const AdaptiveParams& _params,
                              double _dt = config::integrators::default_dt,
                              bool _decent = false)
      : Overdamped(_gamma)
      , Adaptive(_params)
      , IntegratorBase(_dt)
      , decent(_decent)
  {
  }
  ~AdaptiveOverdampedEulerHeun() override = default;

  void step(network& _network)
  {
    auto& nodes = _network.getNodes();
    auto& positions = nodes.positions();
    auto& forces = nodes.forces();

    m_dt = m_nextDt;
    double q = m_params.qMin;
    size_t iter = 0;

    bool error = false;
    bool qGood = false;
    bool energyGood = false;

    double startEnergy = decent ? _network.getEnergy() : 0.0;

    m_rk = positions;
    m_frk = forces;
    while (iter++ < m_params.maxInnerIter) {
      const double overdampedScale = m_dt * m_invGamma;

#pragma omp parallel for schedule(static)
      for (size_t i = 0; i < nodes.size(); i++) {
        positions[i] += forces[i] * overdampedScale;  // eq 8
      }
      // Pos = r_{k+1}bar

      _network.computeForces();
      // forces = f(r_{k+1}bar)

      const double halfOverdampedScale = 0.5 * overdampedScale;
      double estimatedError = -1e10;
#pragma omp parallel for schedule(static) reduction(max : estimatedError)
      for (size_t i = 0; i < nodes.size(); i++) {
        positions[i] = 0.5 * (positions[i] + m_rk[i])
            + forces[i] * halfOverdampedScale;  // eq 9
        // Pos = r_{k+1}

        const double E = (forces[i] - m_frk[i]).norm() * halfOverdampedScale;
        const double tau = m_params.espAbs
            + m_params.espRel * (positions[i] - m_rk[i]).norm();
        estimatedError = std::max(estimatedError, E / tau);
      }
      
      const double estimatedQ = std::pow(0.5 / estimatedError, 2);
      q = std::clamp(estimatedQ, m_params.qMin, m_params.qMax);

      error = (std::isnan(q) || (m_dt == m_params.dtMin && q < 1.0));
      qGood = (q > 1.0);
      energyGood = !decent || (_network.computeEnergy() < startEnergy);

      if ((qGood && energyGood) || error)
        break;

      nodes.positions() = m_rk;
      nodes.forces() = m_frk;
      m_dt = std::clamp(m_dt * q, m_params.dtMin, m_params.dtMax);
    }
    if (error)
      throw std::runtime_error("Adaptive Euler Heun failed to converge");
    m_nextDt = std::clamp(m_dt * q, m_params.dtMin, m_params.dtMax);
  }

private:
  std::vector<Utils::Math::vec2d> m_rk;
  std::vector<Utils::Math::vec2d> m_frk;

  bool decent = false;
};

}  // namespace integration
}  // namespace networkV4
