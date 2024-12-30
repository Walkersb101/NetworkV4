#pragma once

#include <cstdint>

#include <bzlib.h>

namespace networkV4
{

enum class dataOutType : std::uint8_t
{
  CSV,
};

enum class networkOutType : std::uint8_t
{
  None,
  BinV2,
  HDF5,
};

enum class integratorType : std::uint8_t
{
  overdampedEuler,
  overdampedEulerHeun,
  overdampedAdaptiveEulerHeun,
  FireMinimizer,
  OverdampedAdaptiveMinimizer,
};

enum loadVersion : std::uint8_t
{
  BinV1,
  BinV2,
};

enum class StrainType : std::uint8_t
{
  Shear,
  Elongation,
};

enum class protocolType : std::uint8_t
{
  QuasisaticStrain,
  StepStrain,
  Propogator,
};

enum class rootMethod : std::uint8_t
{
  Bisection,
  ITP
};

}  // namespace networkV4