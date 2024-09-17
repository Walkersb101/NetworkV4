#pragma once

#include <array>
#include <type_traits>
#include <stdexcept>
#include <iostream>

// Generic helper to convert enum to its underlying type (integer)
template<typename Enum>
constexpr std::size_t to_index(Enum e) {
    return static_cast<typename std::underlying_type<Enum>::type>(e);
}

// This is a generic map for enum keys
template<typename Enum, typename ValueType, Enum MinEnumValue, Enum MaxEnumValue>
class EnumMap {
public:
    // Determine the number of values in the enum range
    static constexpr std::size_t size = to_index(MaxEnumValue) - to_index(MinEnumValue) + 1;

    // Default constructor initializes the array
    EnumMap() : data_{} {}

    // Accessor for the map (non-const)
    ValueType& operator[](Enum e) {
        check_range(e);
        return data_[to_index(e) - to_index(MinEnumValue)];
    }

    // Accessor for the map (const version)
    const ValueType& operator[](Enum e) const {
        check_range(e);
        return data_[to_index(e) - to_index(MinEnumValue)];
    }

    const auto begin() const { return data_.begin(); }
    const auto end() const { return data_.end(); }

    auto begin() { return data_.begin(); }
    auto end() { return data_.end(); }

private:
    // Static check to ensure the enum value is within the specified range
    static void check_range(Enum e) {
        if (e < MinEnumValue || e > MaxEnumValue) {
            throw std::out_of_range("Enum value out of range");
        }
    }

    // Storage for the map, using an array
    std::array<ValueType, size> data_;
};