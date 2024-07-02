#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>
#include "Bonds.hpp"

TEST_CASE("WLC Bond class", "[bond]")
{
    SECTION("Default constructor")
    {
        networkV4::bond bond;
        REQUIRE(bond.src() == 0);
        REQUIRE(bond.dst() == 0);
        REQUIRE(bond.connected() == false);
        REQUIRE(bond.type() == networkV4::bondType::single);
    }

    SECTION("Constructor with src, dst, connected, type, l0, mu")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        REQUIRE(bond.src() == 1);
        REQUIRE(bond.dst() == 2);
        REQUIRE(bond.connected() == true);
        REQUIRE(bond.type() == networkV4::bondType::sacrificial);
    }

    SECTION("Copy constructor")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        networkV4::bond bondCopy(bond);
        REQUIRE(bondCopy.src() == 1);
        REQUIRE(bondCopy.dst() == 2);
        REQUIRE(bondCopy.connected() == true);
        REQUIRE(bondCopy.type() == networkV4::bondType::sacrificial);
    }

    SECTION("Move constructor")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        networkV4::bond bondMove(std::move(bond));
        REQUIRE(bondMove.src() == 1);
        REQUIRE(bondMove.dst() == 2);
        REQUIRE(bondMove.connected() == true);
        REQUIRE(bondMove.type() == networkV4::bondType::sacrificial);
    }

    SECTION("Copy assignment")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        networkV4::bond bondCopy;
        bondCopy = bond;
        REQUIRE(bondCopy.src() == 1);
        REQUIRE(bondCopy.dst() == 2);
        REQUIRE(bondCopy.connected() == true);
        REQUIRE(bondCopy.type() == networkV4::bondType::sacrificial);
    }

    SECTION("Move assignment")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        networkV4::bond bondMove;
        bondMove = std::move(bond);
        REQUIRE(bondMove.src() == 1);
        REQUIRE(bondMove.dst() == 2);
        REQUIRE(bondMove.connected() == true);
        REQUIRE(bondMove.type() == networkV4::bondType::sacrificial);
    }

    SECTION("Force")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        REQUIRE_THAT(bond.force(1.0), Catch::Matchers::WithinRel(0.0));
    }

    SECTION("Energy")
    {
        networkV4::bond bond(1, 2, true, networkV4::bondType::sacrificial, 1.0, 1.0, 1.0);
        REQUIRE_THAT(bond.energy(1.0), Catch::Matchers::WithinRel(0.0));
    }
}