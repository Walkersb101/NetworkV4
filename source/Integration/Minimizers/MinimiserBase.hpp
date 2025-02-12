#pragma once

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

namespace networkV4
{
namespace minimisation
{

static constexpr double EPS_ENERGY = 1e-8;

struct minimiserParams
{
  double Ftol = config::integrators::miminizer::Ftol;
  double Etol = config::integrators::miminizer::Etol;
  size_t maxIter = config::integrators::miminizer::maxIter;
};

class minimiserBase
{
public:
  minimiserBase() = default;
  minimiserBase(double Ftol, double Etol, size_t maxIter)
      : m_Ftol(Ftol)
      , m_Etol(Etol)
      , m_maxIter(maxIter)
  {
  }
  minimiserBase(const minimiserParams& _params)
      : m_Ftol(_params.Ftol)
      , m_Etol(_params.Etol)
      , m_maxIter(_params.maxIter)
  {
  }
  virtual ~minimiserBase() = default;

public:
  virtual void minimise(network& _network) = 0;

public:
  auto converged(double fdotf, double _Ecurr, double _Eprev) const -> bool
  {
    return (fabs(_Ecurr - _Eprev)
                < m_Etol * 0.5 * (fabs(_Ecurr) + fabs(_Eprev) + EPS_ENERGY)
            || fdotf < m_Ftol * m_Ftol);
  }

public:
  double m_Ftol = config::integrators::miminizer::Ftol;
  double m_Etol = config::integrators::miminizer::Etol;
  size_t m_maxIter = config::integrators::miminizer::maxIter;
};

}  // namespace minimisation
}  // namespace networkV4
