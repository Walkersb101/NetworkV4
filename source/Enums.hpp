#pragma once

#include <cstdint>

#include <bzlib.h>

namespace networkV4
{
enum class bondType : std::uint8_t
{
  single,
  sacrificial,
  matrix,
};

enum class dataOutType : std::uint8_t
{
  CSV,
};

enum class networkOutType : std::uint8_t
{
  None,
  BinV2,
};

enum class intergratorType : std::uint8_t
{
  overdampedEuler,
  overdampedEulerHeun,
  overdampedAdaptiveEulerHeun,
  FireMinimizer,
  OverdampedAdaptiveMinimizer,
};

enum loadVersion : std::uint8_t
{
  binV1 = 1,
  binV2 = 2
};

enum class StrainType : std::uint8_t
{
  Shear,
  Elongation
};

enum class BreakType : std::uint8_t
{
  None,
  Strain,
};

enum class protocolType : std::uint8_t
{
  QuasisaticStrain
};

enum class rootMethod : std::uint8_t
{
  Bisection,
  ITP
};

}  // namespace networkV4