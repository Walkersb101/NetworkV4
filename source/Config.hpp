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
extern double qMin;
extern double qMax;
extern double dtMin;
extern double dtMax;
}  // namespace adaptiveIntergrator

}  // namespace config