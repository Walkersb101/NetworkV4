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

namespace intergrators
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

namespace miminizer
{
extern std::size_t maxIter;
extern double tol;
}  // namespace miminizer

namespace fire2
{
extern double alpha0;
extern std::size_t Ndelay;
extern double finc;
extern double fdec;
extern double falpha;
extern std::size_t Nnegmax;
extern double dmax;
}  // namespace fire2
}  // namespace intergrators

namespace rootMethods
{
extern double targetTol;
extern double minTol;

namespace ITPMethod
{
extern size_t n0;
extern double k1Scale;
extern double k2;
}  // namespace ITPMethod
}  // namespace rootMethods

namespace protocols
{
namespace quasiStaticStrain
{
extern bool errorOnNotSingleBreak;
extern double strainGuessScale;
}  // namespace quasiStaticStrain

namespace stepStrain
{
extern double stressScale;
extern double timeScale;
}  // namespace stepStrain

}  // namespace protocols

namespace IO
{
extern std::string timeDataName;
extern std::string bondDataName;
extern std::string outputType;
extern size_t precision;
}  // namespace IO

namespace hdf5
{
extern std::size_t maxChunk;
extern std::string fileName;
} // namespace hdf5

}  // namespace config