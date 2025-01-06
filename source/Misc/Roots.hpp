#pragma once

#include <cmath>
#include <cstddef>
#include <stdexcept>

#include "Misc/Config.hpp"
#include "Misc/Utils.hpp"

namespace networkV4
{
namespace roots
{
class ITP
{
public:
  ITP() = delete;

  ITP(double _min,
      double _max,
      double _tol,
      size_t _n0 = config::rootMethods::ITPMethod::n0,
      double _k1Scale = config::rootMethods::ITPMethod::k1Scale,
      double _k2 = config::rootMethods::ITPMethod::k2)
      : m_n0(_n0)
      , m_k1Scale(_k1Scale)
      , m_k1(_k1Scale * (_max - _min))
      , m_k2(_k2)
      , m_nHalf(static_cast<size_t>(
            std::ceil(std::log2(abs(_max - _min) / (2 * _tol)))))
      , m_nMax(m_nHalf + _n0)
      , m_tol(_tol)
  {
    if (_min < _max) {
      throw std::invalid_argument("Invalid range");
    }
    if (_tol <= 0.0) {
      throw std::invalid_argument("Invalid tolerance");
    }
  }

  ~ITP() {}

  auto guessRoot(double _a, double _b, double _fa, double _fb) -> double
  {
    double xhalf = (_a + _b) * 0.5;
    double r = m_tol * std::pow(2, m_nMax - m_iters) - (_b - _a) * 0.5;
    double delta = m_k1 * std::pow(_b - _a, m_k2);

    double alpha = 1.0;  // m_iters % 3 == 0 ? 0.5 : 1.0;
    double xf = (alpha * _fb * _a - _fa * _b) / (alpha * _fb - _fa);
    double sigma = Utils::sign(xhalf - xf);
    double xt = delta <= std::abs(xf - xhalf) ? xf + sigma * delta : xhalf;
    double xITP = std::abs(xt - xhalf) <= r ? xt : xhalf - sigma * r;

    m_iters++;
    if (xITP <= _a || xITP >= _b) {
      throw std::invalid_argument("Root not in range");
    }
    return xITP;
  }

  void reset(double _a, double _b)
  {
    m_iters = 0;
    m_nHalf =
        static_cast<size_t>(std::ceil(std::log2(abs(_b - _a) / (2 * m_tol))));
    m_nMax = m_nHalf + m_n0;
  }

  auto nMax() const -> std::size_t { return m_nMax; }

private:
  size_t m_n0;
  double m_k1Scale;
  double m_k1;
  double m_k2;
  std::size_t m_iters = 0;
  std::size_t m_nHalf;
  std::size_t m_nMax;
  double m_tol;
};

}  // namespace roots
}  // namespace networkV4