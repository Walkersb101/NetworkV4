#include "Config.hpp"

namespace config
{
double default_dt = 1e-3;
namespace adaptiveIntergrator
{
double esp = 1e-3;
std::size_t maxIter = 10;

double dtMin = 1e-6;
double dtMax = 1e1;
}  // namespace AdaptiveIntergrator

namespace adaptiveHeun
{
double qMin = 0.001;
double qMax = 1.2;
}  // namespace adaptiveHeun

namespace fire2
{
double alpha0 = 0.1;
std::size_t Ndelay = 5;
double finc = 1.1;
double fdec = 0.5;
double falpha = 0.99;
std::size_t Nnegmax = 2000;
}
}  // namespace config