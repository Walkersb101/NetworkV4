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
  deformBase() = default;
  virtual ~deformBase() = default;

public:
  virtual auto getStrain(const network& _network) const -> double = 0;
  virtual void strain(network& _network, double _step) = 0;
};

class shear : public deformBase
{
public:
  shear() = default;
  virtual ~shear() = default;

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
  elongationAreaY() = default;
  virtual ~elongationAreaY() = default;

public:
  auto getStrain(const network& _network) const -> double override
  {
    return _network.getElongationStrain()[1];
  }

  void strain(network& _network, double _step) override
  {
    double newstrain = _network.getElongationStrain()[1] + _step;

    Utils::Math::vec2d newDomain = Utils::Math::broadcast(
        _network.getRestBox().getDomain(),
        Utils::Math::vec2d({1.0 / (1.0 + newstrain), 1.0 + newstrain}),
                           std::multiplies<>());
    box newBox(newDomain, _network.getBox().getxy());
    _network.setBox(newBox);
  }
};

}  // namespace deform
}  // namespace protocols
}  // namespace networkV4
