#pragma once

#include <string>
#include <unordered_map>
#include <bxzstr.hpp>      // For bxz::Compression
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
using Str2BreakTypeMap = std::unordered_map<std::string, BreakType>;
using BreakType2StrMap = std::unordered_map<BreakType, std::string>;
using Str2ProtocolTypeMap = std::unordered_map<std::string, protocolType>;
using ProtocolType2StrMap = std::unordered_map<protocolType, std::string>;
using Str2RootMethodMap = std::unordered_map<std::string, rootMethod>;
using RootMethod2StrMap = std::unordered_map<rootMethod, std::string>;
using Str2CompressionMap = std::unordered_map<std::string, bxz::Compression>;
using Compression2StrMap = std::unordered_map<bxz::Compression, std::string>;
using CompressionExtMap = std::unordered_map<bxz::Compression, std::string>;

// --- BondType ---
extern const Str2BondTypeMap str2BondType;
extern const BondType2StrMap bondType2Str;

// --- DataOutType ---
extern const Str2DataOutTypeMap str2DataOutType;
extern const DataOutType2StrMap dataOutType2Str;

// --- NetworkOutType ---
extern const Str2NetworkOutTypeMap str2NetworkOutType;
extern const NetworkOutType2StrMap networkOutType2Str;

// --- IntegratorType ---
extern const Str2IntegratorTypeMap str2IntegratorType;
extern const IntegratorType2StrMap integratorType2Str;

// --- LoadVersion ---
extern const Str2LoadVersionMap str2LoadVersion;
extern const LoadVersion2StrMap loadVersion2Str;

// --- StrainType ---
extern const Str2StrainTypeMap str2StrainType;
extern const StrainType2StrMap strainType2Str;

// --- BreakType ---
extern const Str2BreakTypeMap str2BreakType;
extern const BreakType2StrMap breakType2Str;

// --- ProtocolType ---
extern const Str2ProtocolTypeMap str2ProtocolType;
extern const ProtocolType2StrMap protocolType2Str;

// --- RootMethod ---
extern const Str2RootMethodMap str2RootMethod;
extern const RootMethod2StrMap rootMethod2Str;

// --- Compression ---
extern const Str2CompressionMap str2Compression;
extern const Compression2StrMap compression2Str;
extern const CompressionExtMap compressionExt;

}  // namespace enumString
}  // namespace networkV4
