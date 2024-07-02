#pragma once

#include <cstddef>

#define PARALLEL 0

#if PARALLEL
#  include <execution>
#  define PAR std::execution::par,
#else
#  define PAR
#endif

namespace config
{

extern double default_dt;

namespace adaptiveIntergrator
{
extern double esp;
extern std::size_t maxIter;
extern double dtMin;
extern double dtMax;
}  // namespace adaptiveIntergrator

namespace adaptiveHeun
{
extern double qMin;
extern double qMax;
}  // namespace adaptiveHeun

namespace fire2
{
extern double alpha0;
extern size_t Ndelay;
extern double finc;
extern double fdec;
extern double falpha;
extern size_t Nnegmax;
}  // namespace fire2

}  // namespace config