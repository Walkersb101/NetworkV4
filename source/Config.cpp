#include "Config.hpp"

namespace config
{
double default_dt = 1e-3;
namespace adaptiveIntergrator
{
double esp = 1e-3;
std::size_t maxIter = 10;
double qMin = 0.001;
double qMax = 1.2;
double dtMin = 1e-6;
double dtMax = 1e1;
}  // namespace AdaptiveIntergrator
}  // namespace config