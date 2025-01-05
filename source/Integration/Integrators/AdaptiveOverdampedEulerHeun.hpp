#pragma once

#include <algorithm>
#include <vector>

#include <range/v3/view/zip.hpp>

#include "Adaptive.hpp"
#include "Core/Network.hpp"
#include "IntegratorBase.hpp"
#include "Overdamped.hpp"

namespace networkV4
{
namespace integration
{

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
    m_dt = m_nextDt;
    double q = m_params.qMin;
    size_t iter = 0;
    networkV4::nodes& nodes = _network.getNodes();

    bool error = false;
    bool qGood = false;
    bool energyGood = false;

    double startEnergy = decent ? _network.getEnergy() : 0.0;

    m_rn = nodes.positions();
    m_frn = nodes.forces();
    while (iter++ < m_params.maxInnerIter) {
      const double overdampedScale = m_dt * m_invGamma;

      std::transform(nodes.positions().begin(),
                     nodes.positions().end(),
                     nodes.forces().begin(),
                     nodes.positions().begin(),
                     [overdampedScale](const auto& position, const auto& force)
                     { return position + force * overdampedScale; });

      _network.computeForces();
      m_frnbar = nodes.forces();

      const double halfOverdampedScale = 0.5 * overdampedScale;
      for (auto [position, velocity, force, tempForce] : ranges::views::zip(
               nodes.positions(), nodes.velocities(), nodes.forces(), m_frn))
      {
        position += (force - tempForce) * halfOverdampedScale;
        force = (force + tempForce) * 0.5;
        velocity = force * m_invGamma;
      }
      //_network.computeForces();

      const double estimatedError = errorEstimate(nodes);
      const double estimatedQ = std::pow(0.5 / estimatedError, 2);
      q = std::clamp(estimatedQ, m_params.qMin, m_params.qMax);

      error = (std::isnan(q) || (m_dt == m_params.dtMin && q < 1.0));
      q = (q > 1.0);
      energyGood = !decent || (_network.getEnergy() < startEnergy);

      if ((q && energyGood) || error)
        break;

      nodes.positions() = m_rn;
      nodes.forces() = m_frn;
      m_dt = std::clamp(m_dt * q, m_params.dtMin, m_params.dtMax);
    }
    if (error)
      throw std::runtime_error("Adaptive Euler Heun failed to converge");
    m_nextDt = std::clamp(m_dt * q, m_params.dtMin, m_params.dtMax);
    _network.computeForces();
  }

private:
  auto errorEstimate(const nodes& _nodes) const -> double
  {
    // based on https://arxiv.org/pdf/2108.03399
    const double errorScale = m_dt * m_invGamma * 0.5;

    double sum = 0;
    for (auto [frn, frnbar, rn, rnp1] :
         ranges::views::zip(m_frn, m_frnbar, m_rn, _nodes.positions()))
    {
      const double E = (frnbar - frn).norm() * errorScale;
      const double tau = m_params.espAbs + m_params.espRel * (rnp1 - rn).norm();
      const double test = pow(E / tau, 2);
      sum += std::pow(E / tau, 2);
    }
    return std::sqrt(sum);
  }

private:
  std::vector<Utils::vec2d> m_rn;
  std::vector<Utils::vec2d> m_frn;
  std::vector<Utils::vec2d> m_frnbar;

  bool decent = false;
};

}  // namespace integration
}  // namespace networkV4
