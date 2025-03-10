#pragma once

#include <numeric>

#include <range/v3/view/zip.hpp>

#include "Integration/LineSearch/LineSearchQuad.hpp"
#include "MinimiserBase.hpp"

namespace networkV4
{
namespace minimisation

{

class SD : public minimiserBase
{
public:
  SD() = default;
  SD(double Ftol, double Etol, size_t maxIter)
      : minimiserBase(Ftol, Etol, maxIter)
  {
  }

  SD(const minimiserParams& _minParams)
      : minimiserBase(_minParams)
  {
  }

public:
  void minimise(network& _network) override
  {
    auto lineSearch = lineSearch::lineSearchQuad(0.1);

    auto& forces = _network.getNodes().forces();
    _network.computeForces();

    double fdotf = Utils::Math::xdoty(forces, forces);
    if (fdotf < m_Ftol * m_Ftol)
      return;

    double Ecurr = _network.getEnergy();
    double Eprev = Ecurr;

    size_t iter = 0;
    for (iter = 0; iter < m_maxIter; iter++) {
      Eprev = Ecurr;
      auto h = forces;
      auto state = lineSearch.search(h, _network);
      if (!state)
      {  // TODO: If first Step allow more errors
        throw std::runtime_error("Line search failed");  // TODO: Better Logging
      }

      _network.computeForces();
      Ecurr = _network.getEnergy();
      fdotf = Utils::Math::xdoty(forces, forces);
      if (converged(fdotf, Ecurr, Eprev))
        break;
    }
  }
};

}  // namespace minimisation
}  // namespace networkV4