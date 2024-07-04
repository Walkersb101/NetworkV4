#pragma once

#include <cstddef>
#include <cstdint>
#include <string>

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
extern size_t maxIter;
extern double dtMin;
extern double dtMax;
}  // namespace adaptiveIntergrator

namespace adaptiveHeun
{
extern double qMin;
extern double qMax;
}  // namespace adaptiveHeun

namespace miminizer
{
extern size_t maxIter;
extern double tol;
}  // namespace miminizer

namespace fire2
{
extern double alpha0;
extern size_t Ndelay;
extern double finc;
extern double fdec;
extern double falpha;
extern size_t Nnegmax;
extern size_t maxIter;
}  // namespace fire2

namespace ITPMethod
{
extern size_t n0;
extern double k1Scale;
extern double k2;
extern double tol;
}  // namespace ITPMethod

namespace IO{
    extern std::string timeDataName;
    extern std::string bondDataName;
    extern std::string outputType;
}

}  // namespace config