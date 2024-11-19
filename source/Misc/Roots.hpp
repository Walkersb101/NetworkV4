#pragma once
#include <cstddef>

namespace networkV4
{
namespace roots
{
class ITP
{
public:
  ITP();
  ITP(double _min, double _max, double _tol);
  ITP(size_t _n0,
      double _k1Scale,
      double _k2,
      double _min,
      double _max,
      double _tol);
  ~ITP();

public:
  auto guessRoot(double _a, double _b, double _fa, double _fb) -> double;
  void reset(double _a, double _b);

public:
  auto nMax() const -> std::size_t { return m_nMax; }

private:
  size_t m_n0;
  double m_k1Scale;
  double m_k1;
  double m_k2;
  std::size_t m_iters;
  std::size_t m_nHalf;
  std::size_t m_nMax;
  double m_tol;
};

}  // namespace roots
}  // namespace networkV4