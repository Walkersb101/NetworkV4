#pragma once

#include <tl/expected.hpp>

#include "Core/Network.hpp"
#include "LineSearch.hpp"
#include "Misc/Config.hpp"
#include "Misc/Utils.hpp"

namespace quadConfig = config::integrators::lineSearch::quad;

namespace networkV4
{
namespace lineSearch
{

class lineSearchQuad
{
public:
  lineSearchQuad(double _alphaMax)
      : m_alphaMax(_alphaMax)
  {
  }

public:
  void setAlphaMax(double _alphaMax) { m_alphaMax = _alphaMax; }

public:
  auto search(const std::vector<Utils::Math::vec2d>& _h, network& _network)
      -> tl::expected<double, lineSearchState>
  {
    _network.computeForces();
    const double Eoriginal = _network.getEnergy();

    double fdoth = Utils::Math::xdoty(_network.getNodes().forces(), _h);
    if (fdoth <= 0.0)
      return tl::make_unexpected(lineSearchState::DirectionNotDescent);

    double hmax = Utils::maxComp(_h);
    if (hmax < 1e-14)
      return tl::make_unexpected(lineSearchState::zeroforce);

    double alphaMax = std::min(m_alphaMax, quadConfig::alphaMax);

    m_rk = _network.getNodes().positions();
    m_frk = _network.getNodes().forces();

    double alpha = alphaMax;
    double alphaprev = 0.0;
    double fdothprev = fdoth;
    double Ecurr = Eoriginal;
    double Eprev = Ecurr;

    while (true) {
      Ecurr = alphaStep(_network, _h, alpha);

      fdoth = Utils::Math::xdoty(_network.getNodes().forces(), _h);
      double delfh = fdoth - fdothprev;
      if (fabs(fdoth) < quadConfig::esp || fabs(delfh) < quadConfig::esp) {
        resetNetwork(_network);
        return tl::make_unexpected(lineSearchState::zeroquad);
      }

      double relerr = fabs(
          1.0
          - (0.5 * (alpha - alphaprev) * (fdoth + fdothprev) + Ecurr) / Eprev);
      double alpha0 = alpha - (alpha - alphaprev) * fdoth / delfh;
      if (relerr < quadConfig::tol && alpha0 > 0.0 && alpha0 < alphaMax) {
        Ecurr = alphaStep(_network, _h, alpha0);
        if (Ecurr - Eoriginal < quadConfig::EMACH)
          return alpha0;
      }

      double dEIdeal = -quadConfig::backTrackSlope * alpha * fdoth;
      double dE = Ecurr - Eprev;
      if (dE < dEIdeal)
        return alpha;

      fdothprev = fdoth;
      Eprev = Ecurr;
      alphaprev = alpha;

      alpha *= quadConfig::alphaReduce;

      if (alpha <= 0.0 || dEIdeal >= -quadConfig::EMACH) {
        resetNetwork(_network);
        return tl::make_unexpected(lineSearchState::zeroAlpha);
      }
    }
  }

private:
  auto alphaStep(network& _network,
                 const std::vector<Utils::Math::vec2d> _h,
                 const double _alpha) -> double
  {
    std::transform(
        m_rk.begin(),
        m_rk.end(),
        _h.begin(),
        _network.getNodes().positions().begin(),
        [&](const Utils::Math::vec2d& rk, const Utils::Math::vec2d& h)
        { return rk + _alpha * h; });
    _network.computeForces();
    return _network.getEnergy();
  }

  void resetNetwork(network& _network)
  {
    _network.getNodes().positions() = m_rk;
    _network.getNodes().forces() = m_frk;
    _network.computeForces();
  }

private:
  std::vector<Utils::Math::vec2d> m_rk;
  std::vector<Utils::Math::vec2d> m_frk;
  double m_alphaMax = 0.01;
};

}  // namespace lineSearch
}  // namespace networkV4