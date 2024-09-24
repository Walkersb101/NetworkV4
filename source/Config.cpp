#include <cstddef>
#include <cstdint>
#include <string>

#include "Config.hpp"

namespace config
{

namespace intergrators
{
double default_dt = 1e-3;

namespace adaptiveIntergrator
{
double esp = 1e-2;
std::size_t maxIter = 10;

double dtMin = 1e-6;
double dtMax = 5e-1;
}  // namespace adaptiveIntergrator

namespace adaptiveHeun
{
double qMin = 1e-3;
double qMax = 1.2;
}  // namespace adaptiveHeun

namespace miminizer
{
std::size_t maxIter = 1e6;
double tol = 1e-3;
}  // namespace miminizer

namespace fire2
{
double alpha0 = 0.1;
std::size_t Ndelay = 5;
double finc = 1.1;
double fdec = 0.5;
double falpha = 0.99;
std::size_t Nnegmax = 2000;
double dmax = 0.1;
}  // namespace fire2
}  // namespace intergrators

namespace rootMethods
{
double targetTol = 1e-4;
double minTol = 1e-10;
namespace ITPMethod
{
std::size_t n0 = 10;
double k1Scale = 0.2;
double k2 = 2.0;

}  // namespace ITPMethod
}  // namespace rootMethods

namespace protocols
{
namespace quasiStaticStrain
{
bool errorOnNotSingleBreak = false;
double strainGuessScale = 2.0;
}  // namespace quasiStaticStrain

namespace stepStrain
{
double stressScale = 0.9;
double timeScale = 2.0;
}  // namespace stepStrain

}  // namespace protocols

namespace IO
{
std::string timeDataName = "Data";
std::string bondDataName = "Breaks";
std::string outputType = "CSV";
size_t precision = 16;
}  // namespace IO

}  // namespace config