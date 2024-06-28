#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Bonds.hpp"

TEST_CASE("WLC Bond class", "[bond]")
{
    SECTION("Default constructor")
    {
        network::WLCBond bond;
        REQUIRE(bond.src() == 0);
        REQUIRE(bond.dst() == 0);
        REQUIRE(bond.connected() == false);
        REQUIRE(bond.type() == network::bondType::single);
    }

    SECTION("Constructor with src, dst, connected, type, l0, mu")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        REQUIRE(bond.src() == 1);
        REQUIRE(bond.dst() == 2);
        REQUIRE(bond.connected() == true);
        REQUIRE(bond.type() == network::bondType::sacrificial);
    }

    SECTION("Copy constructor")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        network::WLCBond bondCopy(bond);
        REQUIRE(bondCopy.src() == 1);
        REQUIRE(bondCopy.dst() == 2);
        REQUIRE(bondCopy.connected() == true);
        REQUIRE(bondCopy.type() == network::bondType::sacrificial);
    }

    SECTION("Move constructor")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        network::WLCBond bondMove(std::move(bond));
        REQUIRE(bondMove.src() == 1);
        REQUIRE(bondMove.dst() == 2);
        REQUIRE(bondMove.connected() == true);
        REQUIRE(bondMove.type() == network::bondType::sacrificial);
    }

    SECTION("Copy assignment")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        network::WLCBond bondCopy;
        bondCopy = bond;
        REQUIRE(bondCopy.src() == 1);
        REQUIRE(bondCopy.dst() == 2);
        REQUIRE(bondCopy.connected() == true);
        REQUIRE(bondCopy.type() == network::bondType::sacrificial);
    }

    SECTION("Move assignment")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        network::WLCBond bondMove;
        bondMove = std::move(bond);
        REQUIRE(bondMove.src() == 1);
        REQUIRE(bondMove.dst() == 2);
        REQUIRE(bondMove.connected() == true);
        REQUIRE(bondMove.type() == network::bondType::sacrificial);
    }

    SECTION("Force")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        REQUIRE_THAT(bond.force(1.0), Catch::Matchers::WithinRel(0.0));
    }

    SECTION("Energy")
    {
        network::WLCBond bond(1, 2, true, network::bondType::sacrificial, 1.0, 1.0);
        REQUIRE_THAT(bond.energy(1.0), Catch::Matchers::WithinRel(0.0));
    }
}