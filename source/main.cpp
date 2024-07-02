#include "Network.hpp"

auto main() -> int
{
    networkV4::network net;
    net.loadFromBin("/home/sam/Documents/Code/NetworkV4/test/DoubleTest.bin", networkV4::loadVersion::binV2);
    return 0;
}
