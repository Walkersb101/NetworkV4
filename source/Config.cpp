#include "Config.hpp"

#include <cstddef>
#include <cstdint>
#include <string>

namespace config
{

double default_dt = 1e-3;

namespace adaptiveIntergrator
{
double esp = 1e-3;
std::size_t maxIter = 10;

double dtMin = 1e-6;
double dtMax = 5e-1;
}  // namespace adaptiveIntergrator

namespace adaptiveHeun
{
double qMin = 0.001;
double qMax = 1.2;
}  // namespace adaptiveHeun

namespace miminizer
{
std::size_t maxIter = 1e6;
double tol = 1e-4;
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

namespace ITPMethod
{
std::size_t n0 = 20;
double k1Scale = 0.2;
double k2 = 2.0;
double tol = 1e-12;
}  // namespace ITPMethod

namespace IO
{
std::string timeDataName = "Data";
std::string bondDataName = "Breaks";
std::string outputType = "CSV";
size_t precision = 16;
}  // namespace IO

}  // namespace config