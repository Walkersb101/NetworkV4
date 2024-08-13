#include <string>
#include <unordered_map>

#include "EnumString.hpp"

#include <bxzstr.hpp>

#include "Enums.hpp"

namespace networkV4
{
namespace enumString
{
// --- BondType ---
const std::unordered_map<std::string, bondType> str2BondType = {
    {"Single", bondType::single},
    {"Sacrificial", bondType::sacrificial},
    {"Matrix", bondType::matrix}};
const std::unordered_map<bondType, std::string> bondType2Str = {
    {bondType::single, "Single"},
    {bondType::sacrificial, "Sacrificial"},
    {bondType::matrix, "Matrix"}};

// --- DataOutType ---
const std::unordered_map<std::string, dataOutType> str2DataOutType = {
    {"CSV", dataOutType::CSV}};
const std::unordered_map<dataOutType, std::string> dataOutType2Str = {
    {dataOutType::CSV, "CSV"}};

// --- NetworkOutType ---
const std::unordered_map<std::string, networkOutType> str2NetworkOutType = {
    {"None", networkOutType::None}, {"BinV2", networkOutType::BinV2}};
const std::unordered_map<networkOutType, std::string> networkOutType2Str = {
    {networkOutType::None, "None"}, {networkOutType::BinV2, "BinV2"}};

// --- IntergratorType ---
const std::unordered_map<std::string, intergratorType> str2IntergratorType = {
    {"OverdampedEuler", intergratorType::overdampedEuler},
    {"OverdampedEulerHeun", intergratorType::overdampedEulerHeun},
    {"OverdampedAdaptiveEulerHeun",
     intergratorType::overdampedAdaptiveEulerHeun},
    {"FireMinimizer", intergratorType::FireMinimizer},
    {"OverdampedAdaptiveMinimizer",
     intergratorType::OverdampedAdaptiveMinimizer}};
const std::unordered_map<intergratorType, std::string> intergratorType2Str = {
    {intergratorType::overdampedEuler, "OverdampedEuler"},
    {intergratorType::overdampedEulerHeun, "OverdampedEulerHeun"},
    {intergratorType::overdampedAdaptiveEulerHeun,
     "OverdampedAdaptiveEulerHeun"},
    {intergratorType::FireMinimizer, "FireMinimizer"},
    {intergratorType::OverdampedAdaptiveMinimizer,
     "OverdampedAdaptiveMinimizer"}};

// --- LoadVersion ---
const std::unordered_map<std::string, loadVersion> str2LoadVersion = {
    {"BinV1", loadVersion::BinV1}, {"BinV2", loadVersion::BinV2}};
const std::unordered_map<loadVersion, std::string> loadVersion2Str = {
    {loadVersion::BinV1, "BinV1"}, {loadVersion::BinV2, "BinV2"}};

// --- StrainType ---
const std::unordered_map<std::string, StrainType> str2StrainType = {
    {"Shear", StrainType::Shear}, {"Elongation", StrainType::Elongation}};
const std::unordered_map<StrainType, std::string> strainType2Str = {
    {StrainType::Shear, "Shear"}, {StrainType::Elongation, "Elongation"}};

// --- BreakType ---
const std::unordered_map<std::string, BreakType> str2BreakType = {
    {"None", BreakType::None},
    {"Strain", BreakType::Strain},
    {"Energy", BreakType::Energy},
    {"SGR", BreakType::SGR}};
const std::unordered_map<BreakType, std::string> breakType2Str = {
    {BreakType::None, "None"},
    {BreakType::Strain, "Strain"},
    {BreakType::Energy, "Energy"},
    {BreakType::SGR, "SGR"}};

// --- ProtocolType ---
const std::unordered_map<std::string, protocolType> str2ProtocolType = {
    {"QuasiStaticStrain", protocolType::QuasisaticStrain},
    {"StepStrain", protocolType::StepStrain}};
const std::unordered_map<protocolType, std::string> protocolType2Str = {
    {protocolType::QuasisaticStrain, "QuasiStaticStrain"},
    {protocolType::StepStrain, "StepStrain"}};

// --- RootMethod ---
const std::unordered_map<std::string, rootMethod> str2RootMethod = {
    {"Bisection", rootMethod::Bisection}, {"ITP", rootMethod::ITP}};
const std::unordered_map<rootMethod, std::string> rootMethod2Str = {
    {rootMethod::Bisection, "Bisection"}, {rootMethod::ITP, "ITP"}};

// --- Compression ---
const std::unordered_map<std::string, bxz::Compression> str2Compression = {
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
    {"PlainText", bxz::Compression::plaintext}};
const std::unordered_map<bxz::Compression, std::string> compression2Str = {
    {bxz::Compression::z, "Z"},
    {bxz::Compression::bz2, "BZ2"},
    {bxz::Compression::lzma, "LZMA"},
    {bxz::Compression::zstd, "ZSTD"},
    {bxz::Compression::plaintext, "PlainText"}};

const std::unordered_map<bxz::Compression, std::string> compressionExt = {
    {bxz::Compression::z, ".zz"},
    {bxz::Compression::bz2, ".bz2"},
    {bxz::Compression::lzma, ".lzma"},
    {bxz::Compression::zstd, ".zst"},
    {bxz::Compression::plaintext, ""}};

}  // namespace enumString
}  // namespace networkV4
