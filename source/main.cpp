#include "Network.hpp"

auto main() -> int
{
    networkV4::network net;
    net.loadFromBin("network.bin", networkV4::loadVersion::binV2);
    return 0;
}
