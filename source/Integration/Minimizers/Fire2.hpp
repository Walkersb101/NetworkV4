#pragma once

#include <numeric>

#include <range/v3/view/zip.hpp>

#include "MinimiserBase.hpp"
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
              double _dtMax)
      : alpha0(_alpha0)
      , Ndelay(_Ndelay)
      , finc(_finc)
      , fdec(_fdec)
      , falpha(_falpha)
      , Nnegmax(_Nnegmax)
      , dmax(_dmax)
      , dtMin(_dtMin)
      , dtMax(_dtMax)
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

    nodes& nodes = _network.getNodes();
    auto& vels = nodes.velocities();
    auto& forces = nodes.forces();

    _network.computeForces();
    double Eprev = _network.getEnergy();
    double Ecurr = _network.getEnergy();

    // double maxComp = Utils::maxComp(_network.getNodes().forces());

    fdotf = xdoty(forces, forces);
    if (fdotf < m_Ftol * m_Ftol) {
      return;
    }

    nodes.zeroVelocity();
    // updateVelocities(_network, m_dt);

    size_t iter = 0;
    while (iter++ < m_maxIter) {
      vdotf = xdoty(vels, forces);
      if (vdotf > 0.0) {
        Npos++;
        Nneg = 0;

        vdotv = xdoty(vels, vels);
        fdotf = xdoty(forces, forces);

        scale1 = 1.0 - alpha;
        scale2 = fdotf <= 1e-20 ? 0.0 : (alpha * std::sqrt(vdotv / fdotf));

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
        move(_network, -0.5 * m_dt);
        nodes.zeroVelocity();
      }

      updateVelocities(_network, m_dt);

      // double maxAbsVel = maxAbsComponent(vels);
      // m_dt =
      //     std::max(std::min(m_dt, m_params.dmax / maxAbsVel),
      //     m_params.dtMin);

      if (vdotf > 0.0) {
#pragma omp parallel for schedule(static)
        for (size_t i = 0; i < vels.size(); i++) {
          vels[i] *= scale1;
          vels[i] += forces[i] * scale2;
        }
      }
      move(_network, m_dt);

      Eprev = Ecurr;
      _network.computeForces();
      Ecurr = _network.getEnergy();

      if (Npos > m_params.Ndelay
          && fabs(Ecurr - Eprev)
              < m_Etol * 0.5 * (fabs(Ecurr) + fabs(Eprev) + EPS_ENERGY))
      {
        break;
      }

      fdotf = xdoty(forces, forces);
      if (Npos > m_params.Ndelay && fdotf < m_Ftol * m_Ftol) {
        break;
      }
    }
  }

private:
  auto xdoty(const std::vector<Utils::Math::vec2d>& _x,
             const std::vector<Utils::Math::vec2d>& _y) -> double
  {
    double sum = 0.0;
#pragma omp parallel for schedule(static) reduction(+ : sum)
    for (size_t i = 0; i < _x.size(); i++) {
      sum += _x[i] * _y[i];
    }
    return sum;
  }

  auto maxAbsComponent(const std::vector<Utils::Math::vec2d>& _vec) -> double
  {
    return std::accumulate(_vec.begin(),
                           _vec.end(),
                           0.0,
                           [](double _max, const Utils::Math::vec2d& _v)
                           { return std::max(_max, _v.abs().max()); });
  }

  void updateVelocities(network& _network, double _alpha)
  {
    const auto& masses = _network.getNodes().masses();
    auto& vels = _network.getNodes().velocities();
    const auto& forces = _network.getNodes().forces();

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < vels.size(); i++) {
      vels[i] += forces[i] * (_alpha / masses[i]);
    }
  }

  void move(network& _network, double _alpha)
  {
    const auto& vels = _network.getNodes().velocities();
    auto& pos = _network.getNodes().positions();

#pragma omp parallel for schedule(static)
    for (size_t i = 0; i < pos.size(); i++) {
      pos[i] += vels[i] * _alpha;
    }
  }

private:
  Fire2Params m_params = Fire2Params();
  double m_dt = config::integrators::default_dt;
};

}  // namespace minimisation
}  // namespace networkV4