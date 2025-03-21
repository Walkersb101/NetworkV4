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

#define ROUND_ERROR_PRECISION 1.0e-14

#define NUM_TAGS 4

namespace config
{

// Integrators configuration
namespace integrators
{
inline double default_dt = 1e-3;

namespace adaptive
{
inline std::size_t maxIter = 20;

inline double dtMin = 1e-10;
inline double dtMax = 1e-1;

inline double qMin = 1e-3;
inline double qMax = 1.2;

inline double espRel = 0.0;  // 1e-3;
inline double espAbs = 1e-3;

}  // namespace adaptive

namespace lineSearch
{
namespace quad
{
inline double esp = 1e-28;
inline double tol = 0.1;
inline double backTrackSlope = 0.4;
inline double EMACH = 1.0e-8;
inline double alphaMax = 1.0;
inline double alphaReduce = 0.5;
}  // namespace quad
}  // namespace lineSearch

namespace miminizer
{
inline std::size_t maxIter = 1e6;
inline double Ftol = 1e-6;
inline double Etol = 1e-6;
}  // namespace miminizer

namespace fire2
{
inline double alpha0 = 0.25;
inline std::size_t Ndelay = 5;
inline double finc = 1.1;
inline double fdec = 0.5;
inline double falpha = 0.99;
inline std::size_t Nnegmax = 2000;
inline double dmax = 0.1;
}  // namespace fire2

namespace OverdampedAdaptiveMinimizer
{
inline double energyStepScale = 0.5;
inline double dmax = 0.1;
inline double fdec = 0.5;
}  // namespace OverdampedAdaptiveMinimizer
}  // namespace integrators

// Root methods configuration
namespace rootMethods
{
inline double targetTol = 1e-6;
inline double minTol = 1e-10;

namespace ITPMethod
{
inline size_t n0 = 20;
inline double k1Scale = 0.2;
inline double k2 = 2.0;
}  // namespace ITPMethod
}  // namespace rootMethods

// Protocols configuration
namespace protocols
{

inline double maxStep = 1e-2;

namespace quasiStaticStrain
{
inline bool errorOnNotSingleBreak = false;
inline double strainGuessScale = 1.2;
}  // namespace quasiStaticStrain

namespace stepStrain
{
inline double stressScale = 0.9;
inline double timeScale = 2.0;
}  // namespace stepStrain
}  // namespace protocols

// IO configuration
namespace IO
{
namespace timeSeries
{
inline std::string timeDataName = "Data";
inline std::string bondDataName = "Breaks";
}  // namespace timeSeries
namespace CSV
{
inline std::string outputType = "CSV";
inline size_t precision = 16;
}  // namespace CSV

// HDF5 configuration
namespace hdf5
{
inline std::size_t maxChunk = 65536;
inline std::string fileName = "output.h5";
}  // namespace hdf5
}  // namespace IO

namespace partition
{
inline std::size_t mortonRes = 1024;
inline std::size_t passes = 2;
}  // namespace partition

}  // namespace config