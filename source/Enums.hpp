#pragma once

#include <cstdint> // Add this line to include the std::uint8_t type

namespace enums {
    enum class networkType : std::uint8_t { // Change the base type to std::uint8_t
        singleNetwork,
        doubleNetwork
    };
} // namespace enums
