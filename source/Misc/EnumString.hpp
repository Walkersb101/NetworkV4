#pragma once

#include <string>
#include <unordered_map>
#include "Enums.hpp"       // For external enum definitions like bondType, dataOutType, etc.

namespace networkV4 {
namespace enumString {

// --- Type Aliases for Readability ---
using Str2BondTypeMap = std::unordered_map<std::string, bondType>;
using BondType2StrMap = std::unordered_map<bondType, std::string>;
using Str2DataOutTypeMap = std::unordered_map<std::string, dataOutType>;
using DataOutType2StrMap = std::unordered_map<dataOutType, std::string>;
using Str2NetworkOutTypeMap = std::unordered_map<std::string, networkOutType>;
using NetworkOutType2StrMap = std::unordered_map<networkOutType, std::string>;
using Str2IntegratorTypeMap = std::unordered_map<std::string, integratorType>;
using IntegratorType2StrMap = std::unordered_map<integratorType, std::string>;
using Str2LoadVersionMap = std::unordered_map<std::string, loadVersion>;
using LoadVersion2StrMap = std::unordered_map<loadVersion, std::string>;
using Str2StrainTypeMap = std::unordered_map<std::string, StrainType>;
using StrainType2StrMap = std::unordered_map<StrainType, std::string>;
using Str2ProtocolTypeMap = std::unordered_map<std::string, protocolType>;
using ProtocolType2StrMap = std::unordered_map<protocolType, std::string>;
using Str2RootMethodMap = std::unordered_map<std::string, rootMethod>;
using RootMethod2StrMap = std::unordered_map<rootMethod, std::string>;

// --- BondType ---
const Str2BondTypeMap str2BondType = {
    {"Single", bondType::single},
    {"Sacrificial", bondType::sacrificial},
    {"Matrix", bondType::matrix},
    {"Any", bondType::any}
};

const BondType2StrMap bondType2Str = {
    {bondType::single, "Single"},
    {bondType::sacrificial, "Sacrificial"},
    {bondType::matrix, "Matrix"},
    {bondType::any, "Any"}
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
    {"OverdampedEuler", integratorType::overdampedEuler},
    {"OverdampedEulerHeun", integratorType::overdampedEulerHeun},
    {"OverdampedAdaptiveEulerHeun", integratorType::overdampedAdaptiveEulerHeun},
    {"FireMinimizer", integratorType::FireMinimizer},
    {"OverdampedAdaptiveMinimizer", integratorType::OverdampedAdaptiveMinimizer}
};

const IntegratorType2StrMap integratorType2Str = {
    {integratorType::overdampedEuler, "OverdampedEuler"},
    {integratorType::overdampedEulerHeun, "OverdampedEulerHeun"},
    {integratorType::overdampedAdaptiveEulerHeun, "OverdampedAdaptiveEulerHeun"},
    {integratorType::FireMinimizer, "FireMinimizer"},
    {integratorType::OverdampedAdaptiveMinimizer, "OverdampedAdaptiveMinimizer"}
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

// --- ProtocolType ---
const Str2ProtocolTypeMap str2ProtocolType = {
    {"QuasiStaticStrain", protocolType::QuasisaticStrain},
    {"StepStrain", protocolType::StepStrain},
    {"Propogator", protocolType::Propogator}
};

const ProtocolType2StrMap protocolType2Str = {
    {protocolType::QuasisaticStrain, "QuasiStaticStrain"},
    {protocolType::StepStrain, "StepStrain"},
    {protocolType::Propogator, "Propogator"}
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

}  // namespace enumString
}  // namespace networkV4
