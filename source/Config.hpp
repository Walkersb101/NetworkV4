#pragma once

#define PARALLEL 0

#if PARALLEL
#  include <execution>
#  define PAR std::execution::par,
#else
#  define PAR
#endif

struct Config
{
  double default_dt = 1e-3;
};

Config defaultConfig = Config();