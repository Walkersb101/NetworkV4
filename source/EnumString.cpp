#include <string>
#include <unordered_map>

#include "EnumString.hpp"
#include <bxzstr.hpp>
#include "Enums.hpp"

namespace networkV4 {
namespace enumString {

// --- BondType ---
const Str2BondTypeMap str2BondType = {
    {"Single", bondType::single},
    {"Sacrificial", bondType::sacrificial},
    {"Matrix", bondType::matrix}
};

const BondType2StrMap bondType2Str = {
    {bondType::single, "Single"},
    {bondType::sacrificial, "Sacrificial"},
    {bondType::matrix, "Matrix"}
};

// --- DataOutType ---
const Str2DataOutTypeMap str2DataOutType = {
    {"CSV", dataOutType::CSV}
};

const DataOutType2StrMap dataOutType2Str = {
    {dataOutType::CSV, "CSV"}
};

// --- NetworkOutType ---
const Str2NetworkOutTypeMap str2NetworkOutType = {
    {"None", networkOutType::None},
    {"BinV2", networkOutType::BinV2},
    {"HDF5", networkOutType::HDF5}
};

const NetworkOutType2StrMap networkOutType2Str = {
    {networkOutType::None, "None"},
    {networkOutType::BinV2, "BinV2"},
    {networkOutType::HDF5, "HDF5"}
};

// --- IntegratorType ---
const Str2IntegratorTypeMap str2IntegratorType = {
    {"OverdampedEuler", intergratorType::overdampedEuler},
    {"OverdampedEulerHeun", intergratorType::overdampedEulerHeun},
    {"OverdampedAdaptiveEulerHeun", intergratorType::overdampedAdaptiveEulerHeun},
    {"FireMinimizer", intergratorType::FireMinimizer},
    {"OverdampedAdaptiveMinimizer", intergratorType::OverdampedAdaptiveMinimizer}
};

const IntegratorType2StrMap integratorType2Str = {
    {intergratorType::overdampedEuler, "OverdampedEuler"},
    {intergratorType::overdampedEulerHeun, "OverdampedEulerHeun"},
    {intergratorType::overdampedAdaptiveEulerHeun, "OverdampedAdaptiveEulerHeun"},
    {intergratorType::FireMinimizer, "FireMinimizer"},
    {intergratorType::OverdampedAdaptiveMinimizer, "OverdampedAdaptiveMinimizer"}
};

// --- LoadVersion ---
const Str2LoadVersionMap str2LoadVersion = {
    {"BinV1", loadVersion::BinV1},
    {"BinV2", loadVersion::BinV2}
};

const LoadVersion2StrMap loadVersion2Str = {
    {loadVersion::BinV1, "BinV1"},
    {loadVersion::BinV2, "BinV2"}
};

// --- StrainType ---
const Str2StrainTypeMap str2StrainType = {
    {"Shear", StrainType::Shear},
    {"Elongation", StrainType::Elongation}
};

const StrainType2StrMap strainType2Str = {
    {StrainType::Shear, "Shear"},
    {StrainType::Elongation, "Elongation"}
};

// --- BreakType ---
const Str2BreakTypeMap str2BreakType = {
    {"None", BreakType::None},
    {"Strain", BreakType::Strain},
    {"Energy", BreakType::Energy},
    {"SGR", BreakType::SGR}
};

const BreakType2StrMap breakType2Str = {
    {BreakType::None, "None"},
    {BreakType::Strain, "Strain"},
    {BreakType::Energy, "Energy"},
    {BreakType::SGR, "SGR"}
};

// --- ProtocolType ---
const Str2ProtocolTypeMap str2ProtocolType = {
    {"QuasiStaticStrain", protocolType::QuasisaticStrain},
    {"StepStrain", protocolType::StepStrain}
};

const ProtocolType2StrMap protocolType2Str = {
    {protocolType::QuasisaticStrain, "QuasiStaticStrain"},
    {protocolType::StepStrain, "StepStrain"}
};

// --- RootMethod ---
const Str2RootMethodMap str2RootMethod = {
    {"Bisection", rootMethod::Bisection},
    {"ITP", rootMethod::ITP}
};

const RootMethod2StrMap rootMethod2Str = {
    {rootMethod::Bisection, "Bisection"},
    {rootMethod::ITP, "ITP"}
};

// --- Compression ---
const Str2CompressionMap str2Compression = {
#if BXZSTR_Z_SUPPORT
    {"Z", bxz::Compression::z},
#endif
#if BXZSTR_BZ2_SUPPORT
    {"BZ2", bxz::Compression::bz2},
#endif
#if BXZSTR_LZMA_SUPPORT
    {"LZMA", bxz::Compression::lzma},
#endif
#if BXZSTR_ZSTD_SUPPORT
    {"ZSTD", bxz::Compression::zstd},
#endif
    {"PlainText", bxz::Compression::plaintext}
};

const Compression2StrMap compression2Str = {
    {bxz::Compression::z, "Z"},
    {bxz::Compression::bz2, "BZ2"},
    {bxz::Compression::lzma, "LZMA"},
    {bxz::Compression::zstd, "ZSTD"},
    {bxz::Compression::plaintext, "PlainText"}
};

const CompressionExtMap compressionExt = {
    {bxz::Compression::z, ".zz"},
    {bxz::Compression::bz2, ".bz2"},
    {bxz::Compression::lzma, ".lzma"},
    {bxz::Compression::zstd, ".zst"},
    {bxz::Compression::plaintext, ""}
};

}  // namespace enumString
}  // namespace networkV4
