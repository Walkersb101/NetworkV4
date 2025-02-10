#pragma once

#include <numeric>

#include <range/v3/view/zip.hpp>

#include "MinimiserBase.hpp"
#include "Misc/Math/Misc.hpp"
#include "Misc/Utils.hpp"

namespace networkV4
{
namespace minimisation
{

struct Fire2Params
{
  Fire2Params() = default;
  Fire2Params(double _alpha0,
              size_t _Ndelay,
              double _finc,
              double _fdec,
              double _falpha,
              size_t _Nnegmax,
              double _dmax,
              double _dtMin,
              double _dtMax,
              bool _abc)
      : alpha0(_alpha0)
      , Ndelay(_Ndelay)
      , finc(_finc)
      , fdec(_fdec)
      , falpha(_falpha)
      , Nnegmax(_Nnegmax)
      , dmax(_dmax)
      , dtMin(_dtMin)
      , dtMax(_dtMax)
      , abc(_abc)
  {
  }
  double alpha0 = config::integrators::fire2::alpha0;
  size_t Ndelay = config::integrators::fire2::Ndelay;
  double finc = config::integrators::fire2::finc;
  double fdec = config::integrators::fire2::fdec;
  double falpha = config::integrators::fire2::falpha;
  size_t Nnegmax = config::integrators::fire2::Nnegmax;
  double dmax = config::integrators::fire2::dmax;

  double dtMin = config::integrators::adaptive::dtMin;
  double dtMax = config::integrators::adaptive::dtMax;

  bool abc = true;
};

class fire2 : public minimiserBase
{
public:
  fire2() = default;
  fire2(double Ftol,
        double Etol,
        size_t maxIter,
        const Fire2Params& _params,
        double _dt = config::integrators::default_dt)
      : minimiserBase(Ftol, Etol, maxIter)
      , m_params(_params)
      , m_dt(_dt)
  {
  }
  fire2(const minimiserParams& _minParams,
        const Fire2Params& _params = Fire2Params(),
        double _dt = config::integrators::default_dt)
      : minimiserBase(_minParams)
      , m_params(_params)
      , m_dt(_dt)
  {
  }

public:
  void minimise(network& _network) override
  {
    size_t Npos = 0;
    size_t Nneg = 0;
    double alpha = m_params.alpha0;

    double scale1, scale2, vdotf, vdotv, fdotf;

    const auto& masses = _network.getNodes().masses();
    const auto& forces = _network.getNodes().forces();

    nodes& nodes = _network.getNodes();
    auto& vels = nodes.velocities();
    auto& pos = _network.getNodes().positions();

    _network.computeForces();
    double Eprev = _network.getEnergy();
    double Ecurr = _network.getEnergy();

    fdotf = Utils::Math::xdoty(forces, forces);
    if (fdotf < m_Ftol * m_Ftol) {
      return;
    }

    nodes.zeroVelocity();

    size_t iter = 0;
    while (iter++ < m_maxIter) {
      vdotf = Utils::Math::xdoty(vels, forces);

      if (vdotf > 0.0) {
        Npos++;
        Nneg = 0;

        vdotv = Utils::Math::xdoty(vels, vels);
        fdotf = Utils::Math::xdoty(forces, forces);

        if (m_params.abc) {
          alpha = std::max(alpha, 1e-10);
          double abc = 1.0 - std::pow(1.0 - alpha, Npos);
          scale1 = (1.0 - alpha) / abc;
          scale2 =
              fdotf <= 1e-20 ? 0.0 : (alpha * std::sqrt(vdotv / fdotf)) / abc;
        } else {
          scale1 = 1.0 - alpha;
          scale2 = fdotf <= 1e-20 ? 0.0 : (alpha * std::sqrt(vdotv / fdotf));
        }

        if (Npos > m_params.Ndelay) {
          m_dt = std::min(m_dt * m_params.finc, m_params.dtMax);
          alpha *= m_params.falpha;
        }
      } else {
        Nneg++;
        Npos = 0;

        if (Nneg > m_params.Nnegmax) {
          break;
        }
        if (iter > m_params.Ndelay) {
          m_dt = std::max(m_dt * m_params.fdec, m_params.dtMin);
          alpha = m_params.alpha0;
        }

        double scale = 0.5 * m_dt;
        std::transform(pos.begin(),
                       pos.end(),
                       vels.begin(),
                       pos.begin(),
                       [scale](auto& p, auto& v) { return p - scale * v; });
        nodes.zeroVelocity();
      }

      for (size_t i = 0; i < vels.size(); i++) {
        vels[i] += (m_dt / masses[i]) * forces[i];
      }
      if (vdotf > 0.0) {
        std::transform(vels.begin(),
                       vels.end(),
                       forces.begin(),
                       vels.begin(),
                       [scale1, scale2](auto& v, auto& f)
                       { return scale1 * v + scale2 * f; });
        if (m_params.abc) {
          double vmax = m_params.dmax / m_dt;
          std::transform(vels.begin(),
                        vels.end(),
                        vels.begin(),
                        [vmax](auto& v)
                        { return Utils::clamp(v, -vmax, vmax); });

        }
      }
      std::transform(pos.begin(),
                     pos.end(),
                     vels.begin(),
                     pos.begin(),
                     [this](auto& p, auto& v) { return p + m_dt * v; });

      Eprev = Ecurr;
      _network.computeForces();
      Ecurr = _network.getEnergy();

      if (Npos > m_params.Ndelay
          && fabs(Ecurr - Eprev)
              < m_Etol * 0.5 * (fabs(Ecurr) + fabs(Eprev) + EPS_ENERGY))
      {
        break;
      }

      fdotf = Utils::Math::xdoty(forces, forces);
      if (Npos > m_params.Ndelay && fdotf < m_Ftol * m_Ftol) {
        break;
      }
    }
  }

  // auto maxAbsComponent(const std::vector<Utils::Math::vec2d>& _vec) -> double
  //{
  //   return std::accumulate(_vec.begin(),
  //                          _vec.end(),
  //                          0.0,
  //                          [](double _max, const Utils::Math::vec2d& _v)
  //                          { return std::max(_max, _v.abs().max()); });
  // }

private:
  Fire2Params m_params = Fire2Params();
  double m_dt = config::integrators::default_dt;
};

}  // namespace minimisation
}  // namespace networkV4