#pragma once

#include <numeric>

#include <range/v3/view/zip.hpp>

#include "Integration/Integrators/AdaptiveOverdampedEulerHeun.hpp"
#include "Integration/RunStyles/Run.hpp"
#include "MinimiserBase.hpp"

namespace networkV4
{
namespace minimisation
{

class AdaptiveHeunDecent : public minimiserBase
{
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
    integration::AdaptiveOverdampedEulerHeun stepper(1.0, m_params, m_dt, true);

    _network.computeForces();

    double Ecurr = _network.getEnergy();
    double Eprev = Ecurr;

    size_t iter = 0;
    while (iter++ < m_maxIter) {
      Eprev = Ecurr;
      stepper.step(_network);
      _network.computeForces();
      Ecurr = _network.getEnergy();

      if (fabs(Ecurr - Eprev)
          < m_Etol * 0.5 * (fabs(Ecurr) + fabs(Eprev) + EPS_ENERGY))
        return;

      double fdotf =
          xdoty(_network.getNodes().forces(), _network.getNodes().forces());
      if (fdotf < m_Ftol * m_Ftol) {
        return;
      }
    }
  };

private:
  auto xdoty(const std::vector<Utils::vec2d>& _x,
             const std::vector<Utils::vec2d>& _y) -> double
  {
    return std::inner_product(_x.begin(),
                              _x.end(),
                              _y.begin(),
                              0.0,
                              std::plus<double>(),
                              [](Utils::vec2d _a, Utils::vec2d _b)
                              { return _a.dot(_b); });
  }

private:
  integration::AdaptiveParams m_params = integration::AdaptiveParams();
  double m_dt = config::integrators::default_dt;
};

}  // namespace minimisation
}  // namespace networkV4