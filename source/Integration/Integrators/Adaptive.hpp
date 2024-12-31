#pragma once

#pragma once

#include <algorithm>
#include <vector>

#include "Misc/Config.hpp"

namespace networkV4
{
namespace integration
{

struct AdaptiveParams
{
  AdaptiveParams(size_t _maxInnerIter,
                 double _dtMin,
                 double _dtMax,
                 double _qMin,
                 double _qMax,
                 double _espRel,
                 double _espAbs)
      : maxInnerIter(_maxInnerIter)
      , dtMin(_dtMin)
      , dtMax(_dtMax)
      , qMin(_qMin)
      , qMax(_qMax)
      , espRel(_espRel)
      , espAbs(_espAbs)
  {
  }
  size_t maxInnerIter = config::integrators::adaptive::maxIter;
  double dtMin = config::integrators::adaptive::dtMin;
  double dtMax = config::integrators::adaptive::dtMax;
  double qMin = config::integrators::adaptive::qMin;
  double qMax = config::integrators::adaptive::qMax;
  double espRel = config::integrators::adaptive::espRel;
  double espAbs = config::integrators::adaptive::espAbs;
};

class Adaptive
{
public:
  Adaptive();
  Adaptive(AdaptiveParams _params)
      : m_params(_params)
  {
  }
  virtual ~Adaptive() = default;

public:
    auto getNexDt() const -> double { return m_nextDt; }

protected:
  AdaptiveParams m_params;

  double m_nextDt = config::integrators::default_dt;
};

}  // namespace integration
}  // namespace networkV4