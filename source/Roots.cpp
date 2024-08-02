#include <cassert>
#include <cmath>

#include "Roots.hpp"

#include "Config.hpp"
#include "Tools.hpp"

networkV4::roots::ITP::ITP()
    : m_n0(config::rootMethods::ITPMethod::n0)
    , m_k1Scale(config::rootMethods::ITPMethod::k1Scale)
    , m_k2(config::rootMethods::ITPMethod::k2)
    , m_iters(0)
    , m_nHalf(0)
    , m_nMax(0)
    , m_tol(config::rootMethods::targetTol)
{
}

networkV4::roots::ITP::ITP(double _min, double _max, double _tol)
    : ITP(config::rootMethods::ITPMethod::n0,
          config::rootMethods::ITPMethod::k1Scale,
          config::rootMethods::ITPMethod::k2,
          _min,
          _max,
          _tol)
{
}

networkV4::roots::ITP::ITP(size_t _n0,
                           double _k1Scale,
                           double _k2,
                           double _min,
                           double _max,
                           double _tol)
    : m_n0(_n0)
    , m_k1Scale(_k1Scale)
    , m_k1(_k1Scale * (_max - _min))
    , m_k2(_k2)
    , m_iters(0)
    , m_nHalf(static_cast<size_t>(std::ceil(std::log2((_max - _min) / (2 * _tol)))))
    , m_nMax(m_nHalf + _n0)
    , m_tol(_tol)
{
  if (_min < _max) {
    std::invalid_argument("Invalid range");
  }
  if (_tol <= 0.0) {
    std::invalid_argument("Invalid tolerance");
  }
}

networkV4::roots::ITP::~ITP() {}

auto networkV4::roots::ITP::guessRoot(double _a,
                                      double _b,
                                      double _fa,
                                      double _fb) -> double
{
  double xhalf = (_a + _b) * 0.5;
  double r = m_tol * std::pow(2, m_nMax - m_iters) - (_b - _a) * 0.5;
  double delta = m_k1 * std::pow(_b - _a, m_k2);

  double alpha = 1.0;//m_iters % 3 == 0 ? 0.5 : 1.0;
  double xf = (alpha * _fb * _a - _fa * _b) / (alpha * _fb - _fa);
  double sigma = tools::sign(xhalf - xf);
  double xt = delta <= std::abs(xf - xhalf) ? xf + sigma * delta : xhalf;
  double xITP = std::abs(xt - xhalf) <= r ? xt : xhalf - sigma * r;

  m_iters++;
  if (xITP <= _a || xITP >= _b ){
    std::invalid_argument("Root not in range");
  }
  return xITP;
}

void networkV4::roots::ITP::reset(double _min, double _max)
{
  assert(_min < _max);
  m_iters = 0;
  m_nHalf = static_cast<size_t>(std::ceil(std::log2((_max - _min) / (2 * m_tol))));
  m_nMax = m_nHalf + m_n0;
}