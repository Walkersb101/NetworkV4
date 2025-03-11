#pragma once

#include "Partition.hpp"

namespace networkV4
{
namespace OMP
{

extern partition::Partitions threadPartitions;
extern size_t passes;
extern std::vector<networkV4::stresses> localStresses;
extern std::vector<networkV4::bondQueue> localBreaks;

} // namespace OMP
}  // namespace networkV4