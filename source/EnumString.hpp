#pragma once

#include <string>
#include <unordered_map>

#include <bxzstr.hpp>

#include "Enums.hpp"

namespace networkV4
{
namespace enumString
{
// --- BondType ---
extern const std::unordered_map<std::string, bondType> str2BondType;
extern const std::unordered_map<bondType, std::string> bondType2Str;

// --- DataOutType ---
extern const std::unordered_map<std::string, dataOutType> str2DataOutType;
extern const std::unordered_map<dataOutType, std::string> dataOutType2Str;

// --- NetworkOutType ---
extern const std::unordered_map<std::string, networkOutType> str2NetworkOutType;
extern const std::unordered_map<networkOutType, std::string> networkOutType2Str;

// --- IntergratorType ---
extern const std::unordered_map<std::string, intergratorType>
    str2IntergratorType;
extern const std::unordered_map<intergratorType, std::string>
    intergratorType2Str;

// --- LoadVersion ---
extern const std::unordered_map<std::string, loadVersion> str2LoadVersion;
extern const std::unordered_map<loadVersion, std::string> loadVersion2Str;

// --- StrainType ---
extern const std::unordered_map<std::string, StrainType> str2StrainType;
extern const std::unordered_map<StrainType, std::string> strainType2Str;

// --- BreakType ---
extern const std::unordered_map<std::string, BreakType> str2BreakType;
extern const std::unordered_map<BreakType, std::string> breakType2Str;

// --- ProtocolType ---
extern const std::unordered_map<std::string, protocolType> str2ProtocolType;
extern const std::unordered_map<protocolType, std::string> protocolType2Str;

// --- RootMethod ---
extern const std::unordered_map<std::string, rootMethod> str2RootMethod;
extern const std::unordered_map<rootMethod, std::string> rootMethod2Str;

// --- Compression ---
extern const std::unordered_map<std::string, bxz::Compression> str2Compression;
extern const std::unordered_map<bxz::Compression, std::string> compression2Str;
extern const std::unordered_map<bxz::Compression, std::string> compressionExt;
}  // namespace enumString
}  // namespace networkV4
