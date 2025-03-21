#pragma once

#include <numeric>

#include <range/v3/view/zip.hpp>

#include "Integration/Integrators/Overdamped/AdaptiveEulerHeun.hpp"
#include "Integration/RunStyles/Run.hpp"
#include "MinimiserBase.hpp"

namespace networkV4
{
namespace minimisation
{

class AdaptiveHeunDecent : public minimiserBase
{
public:
  AdaptiveHeunDecent() = default;
  AdaptiveHeunDecent(double Ftol,
                     double Etol,
                     size_t maxIter,
                     const integration::AdaptiveParams& _params,
                     double _h = config::integrators::default_dt)
      : minimiserBase(Ftol, Etol, maxIter)
      , m_params(_params)
      , m_dt(_h)
  {
  }

  AdaptiveHeunDecent(const minimiserParams& _minParams,
                     const integration::AdaptiveParams& _params,
                     double _dt = config::integrators::default_dt)
      : minimiserBase(_minParams)
      , m_params(_params)
      , m_dt(_dt)
  {
  }

public:
  void minimise(network& _network) override
  {
    integration::AdaptiveOverdampedEulerHeun stepper(1.0, m_params, m_dt);

    _network.computeForces();

    double Ecurr = _network.getEnergy();
    double Eprev = Ecurr;

    double fdotf = Utils::Math::xdoty(_network.getNodes().forces(),
                                      _network.getNodes().forces());
    if (fdotf < m_Ftol * m_Ftol)
      return;

    double ediff = 0.0;

    size_t iter = 0;
    while (iter++ < m_maxIter) {
      Eprev = Ecurr;
      stepper.step(_network);
      _network.computeForces();
      Ecurr = _network.getEnergy();

      fdotf = Utils::Math::xdoty(_network.getNodes().forces(),
                                 _network.getNodes().forces());
      if (converged(fdotf, Ecurr, Eprev))
        return;
    }
  };

private:
  integration::AdaptiveParams m_params = integration::AdaptiveParams();
  double m_dt = config::integrators::default_dt;
};

}  // namespace minimisation
}  // namespace networkV4