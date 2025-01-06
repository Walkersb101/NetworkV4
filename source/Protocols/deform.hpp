#pragma once

#include "Core/Network.hpp"

namespace networkV4
{
namespace protocols
{
namespace deform
{

class deformBase
{
public:
  deformBase();
  virtual ~deformBase();

public:
  virtual auto getStrain(const network& _network) const -> double = 0;
  virtual void strain(network& _network, double _step) = 0;
};

class shear : public deformBase
{
public:
  shear();
  virtual ~shear();

public:
  auto getStrain(const network& _network) const -> double override
  {
    return _network.getShearStrain();
  }

  void strain(network& _network, double _step) override
  {
    _network.shear(_step);
  }
};

class elongationAreaY : public deformBase
{
public:
  elongationAreaY();
  virtual ~elongationAreaY();

public:
  auto getStrain(const network& _network) const -> double override
  {
    return _network.getElongationStrain().y;
  }

  void strain(network& _network, double _step) override
  {
    double newstrain = _network.getElongationStrain().y + _step;

    Utils::vec2d newDomain = _network.getRestBox().getDomain()
        * Utils::vec2d(1.0 / (1.0 + newstrain), 1.0 + newstrain);
    box newBox(newDomain, _network.getBox().getxy());
    _network.setBox(newBox);
  }
};

}  // namespace deform
}  // namespace protocols
}  // namespace networkV4
