#pragma once

#if PARALLEL
#include <execution>
#define PAR std::execution::par,
#else
#define PAR
#endif

#define OUTPRECISION 16