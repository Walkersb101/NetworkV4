#pragma once

#include "Core/Network.hpp"
#include "Misc/Config.hpp"

namespace networkV4
{
namespace minimisation
{

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
  virtual ~minimiserBase() = default;

public:
  virtual void minimise(network& _network) = 0;

public:
  double m_Ftol = config::integrators::miminizer::Ftol;
  double m_Etol = config::integrators::miminizer::Etol;
  size_t m_maxIter = config::integrators::miminizer::maxIter;
};

}  // namespace integration
}  // namespace networkV4
