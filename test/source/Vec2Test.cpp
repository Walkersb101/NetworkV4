#include "Vec2.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("Vec2d Constructors", "[Vec2]")
{
  SECTION("Default Constructor")
  {
    vec2<double> v;
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(0.0));
  }

  SECTION("Constructor")
  {
    vec2<double> v(2.0, 3.0);
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(2.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(3.0));
  }

  SECTION("Copy Constructor")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(v1);
    REQUIRE_THAT(v2.x, Catch::Matchers::WithinRel(2.0));
    REQUIRE_THAT(v2.y, Catch::Matchers::WithinRel(3.0));
  }

  SECTION("Assignment Operator")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2;
    v2 = v1;
    REQUIRE_THAT(v2.x, Catch::Matchers::WithinRel(2.0));
    REQUIRE_THAT(v2.y, Catch::Matchers::WithinRel(3.0));
  }
}

TEST_CASE("Vec2d Scalar Operators", "[Vec2]")
{
  SECTION("Scalar Addition")
  {
    vec2<double> v(2.0, 3.0);
    vec2<double> result = v + 2.0;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(5.0));
  }

  SECTION("Scalar Subtraction")
  {
    vec2<double> v(2.0, 3.0);
    vec2<double> result = v - 2.0;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(0.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(1.0));
  }

  SECTION("Scalar Multiplication")
  {
    vec2<double> v(2.0, 3.0);
    vec2<double> result = v * 2.0;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(6.0));
  }

  SECTION("Scalar Division")
  {
    vec2<double> v(8.0, 12.0);
    vec2<double> result = v / 2.0;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(6.0));
  }
}

TEST_CASE("Vec2d Vector Operators", "[Vec2]")
{
  SECTION("Addition")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    vec2<double> result = v1 + v2;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(6.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(8.0));
  }

  SECTION("Subtraction")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    vec2<double> result = v1 - v2;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(-2.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(-2.0));
  }

  SECTION("Multiplication")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    vec2<double> result = v1 * v2;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(8.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(15.0));
  }

  SECTION("Division")
  {
    vec2<double> v1(8.0, 12.0);
    vec2<double> v2(2.0, 3.0);
    vec2<double> result = v1 / v2;
    REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(4.0));
  }
}

TEST_CASE("Vec2d Scalar Operator Assigment", "[Vec2]")
{
  SECTION("Scalar Addition Assignment")
  {
    vec2<double> v(2.0, 3.0);
    v += 2.0;
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(5.0));
  }

  SECTION("Scalar Subtraction Assignment")
  {
    vec2<double> v(2.0, 3.0);
    v -= 2.0;
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(1.0));
  }

  SECTION("Scalar Multiplication Assignment")
  {
    vec2<double> v(2.0, 3.0);
    v *= 2.0;
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(6.0));
  }

  SECTION("Scalar Division Assignment")
  {
    vec2<double> v(8.0, 12.0);
    v /= 2.0;
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(6.0));
  }
}

TEST_CASE("Vec2d Vector Operator Assigment", "[Vec2]")
{
  SECTION("Addition Assignment")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    v1 += v2;
    REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(6.0));
    REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(8.0));
  }

  SECTION("Subtraction Assignment")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    v1 -= v2;
    REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(-2.0));
    REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(-2.0));
  }

  SECTION("Multiplication Assignment")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    v1 *= v2;
    REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(8.0));
    REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(15.0));
  }

  SECTION("Division Assignment")
  {
    vec2<double> v1(8.0, 12.0);
    vec2<double> v2(2.0, 3.0);
    v1 /= v2;
    REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(4.0));
  }
}

TEST_CASE("Vec2d Set", "[Vec2]")
{
  SECTION("Set")
  {
    vec2<double> v;
    v.set(2.0, 3.0);
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(2.0));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(3.0));
  }

}

TEST_CASE("Vec2d Lin Alg", "[Vec2]")
{
  SECTION("Normalize")
  {
    vec2<double> v(3.0, 4.0);
    v.normalize();
    REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.6));
    REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(0.8));
  }

  SECTION("Length")
  {
    vec2<double> v(3.0, 4.0);
    double length = v.length();
    REQUIRE_THAT(length, Catch::Matchers::WithinRel(5.0));
  }

  SECTION("Length Squared")
  {
    vec2<double> v(3.0, 4.0);
    double lengthSquared = v.lengthSquared();
    REQUIRE_THAT(lengthSquared, Catch::Matchers::WithinRel(25.0));
  }

  SECTION("Dot")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    double dot = v1.dot(v2);
    REQUIRE_THAT(dot, Catch::Matchers::WithinRel(23.0));
  }

  SECTION("Cross")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    double cross = v1.cross(v2);
    REQUIRE_THAT(cross, Catch::Matchers::WithinRel(-2.0));
  }
}

TEST_CASE("Vec2d Min Max", "[Vec2]")
{
  SECTION("Max")
  {
    vec2<double> v(2.0, 3.0);
    double max = v.max();
    REQUIRE_THAT(max, Catch::Matchers::WithinRel(3.0));
  }

  SECTION("Min")
  {
    vec2<double> v(2.0, 3.0);
    double min = v.min();
    REQUIRE_THAT(min, Catch::Matchers::WithinRel(2.0));
  }

  SECTION("Max Vec2")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    vec2<double> max = v1.max(v2);
    REQUIRE_THAT(max.x, Catch::Matchers::WithinRel(4.0));
    REQUIRE_THAT(max.y, Catch::Matchers::WithinRel(5.0));
  }

  SECTION("Min Vec2")
  {
    vec2<double> v1(2.0, 3.0);
    vec2<double> v2(4.0, 5.0);
    vec2<double> min = v1.min(v2);
    REQUIRE_THAT(min.x, Catch::Matchers::WithinRel(2.0));
    REQUIRE_THAT(min.y, Catch::Matchers::WithinRel(3.0));
  }
}

TEST_CASE("Vec2f Constructors", "[Vec2]")
{
    SECTION("Default Constructor")
    {
        vec2<float> v;
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(0.0f));
    }

    SECTION("Constructor")
    {
        vec2<float> v(2.0f, 3.0f);
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(2.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(3.0f));
    }

    SECTION("Copy Constructor")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(v1);
        REQUIRE_THAT(v2.x, Catch::Matchers::WithinRel(2.0f));
        REQUIRE_THAT(v2.y, Catch::Matchers::WithinRel(3.0f));
    }

    SECTION("Assignment Operator")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2;
        v2 = v1;
        REQUIRE_THAT(v2.x, Catch::Matchers::WithinRel(2.0f));
        REQUIRE_THAT(v2.y, Catch::Matchers::WithinRel(3.0f));
    }
}
TEST_CASE("Vec2f Scalar Operators", "[Vec2]")
{
    SECTION("Scalar Addition")
    {
        vec2<float> v(2.0f, 3.0f);
        vec2<float> result = v + 2.0f;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(5.0f));
    }

    SECTION("Scalar Subtraction")
    {
        vec2<float> v(2.0f, 3.0f);
        vec2<float> result = v - 2.0f;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(0.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(1.0f));
    }

    SECTION("Scalar Multiplication")
    {
        vec2<float> v(2.0f, 3.0f);
        vec2<float> result = v * 2.0f;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(6.0f));
    }

    SECTION("Scalar Division")
    {
        vec2<float> v(8.0f, 12.0f);
        vec2<float> result = v / 2.0f;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(6.0f));
    }
}
TEST_CASE("Vec2f Vector Operators", "[Vec2]")
{
    SECTION("Addition")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        vec2<float> result = v1 + v2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(6.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(8.0f));
    }

    SECTION("Subtraction")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        vec2<float> result = v1 - v2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(-2.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(-2.0f));
    }

    SECTION("Multiplication")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        vec2<float> result = v1 * v2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(8.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(15.0f));
    }

    SECTION("Division")
    {
        vec2<float> v1(8.0f, 12.0f);
        vec2<float> v2(2.0f, 3.0f);
        vec2<float> result = v1 / v2;
        REQUIRE_THAT(result.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(result.y, Catch::Matchers::WithinRel(4.0f));
    }
}
TEST_CASE("Vec2f Scalar Operator Assigment", "[Vec2]")
{
    SECTION("Scalar Addition Assignment")
    {
        vec2<float> v(2.0f, 3.0f);
        v += 2.0f;
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(5.0f));
    }

    SECTION("Scalar Subtraction Assignment")
    {
        vec2<float> v(2.0f, 3.0f);
        v -= 2.0f;
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(1.0f));
    }

    SECTION("Scalar Multiplication Assignment")
    {
        vec2<float> v(2.0f, 3.0f);
        v *= 2.0f;
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(6.0f));
    }

    SECTION("Scalar Division Assignment")
    {
        vec2<float> v(8.0f, 12.0f);
        v /= 2.0f;
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(6.0f));
    }
}
TEST_CASE("Vec2f Vector Operator Assigment", "[Vec2]")
{
    SECTION("Addition Assignment")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        v1 += v2;
        REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(6.0f));
        REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(8.0f));
    }

    SECTION("Subtraction Assignment")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        v1 -= v2;
        REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(-2.0f));
        REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(-2.0f));
    }

    SECTION("Multiplication Assignment")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        v1 *= v2;
        REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(8.0f));
        REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(15.0f));
    }

    SECTION("Division Assignment")
    {
        vec2<float> v1(8.0f, 12.0f);
        vec2<float> v2(2.0f, 3.0f);
        v1 /= v2;
        REQUIRE_THAT(v1.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(v1.y, Catch::Matchers::WithinRel(4.0f));
    }
}
TEST_CASE("Vec2f Set", "[Vec2]")
{
    SECTION("Set")
    {
        vec2<float> v;
        v.set(2.0f, 3.0f);
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(2.0f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(3.0f));
    }

}
TEST_CASE("Vec2f Lin Alg", "[Vec2]")
{
    SECTION("Normalize")
    {
        vec2<float> v(3.0f, 4.0f);
        v.normalize();
        REQUIRE_THAT(v.x, Catch::Matchers::WithinRel(0.6f));
        REQUIRE_THAT(v.y, Catch::Matchers::WithinRel(0.8f));
    }

    SECTION("Length")
    {
        vec2<float> v(3.0f, 4.0f);
        float length = v.length();
        REQUIRE_THAT(length, Catch::Matchers::WithinRel(5.0f));
    }

    SECTION("Length Squared")
    {
        vec2<float> v(3.0f, 4.0f);
        float lengthSquared = v.lengthSquared();
        REQUIRE_THAT(lengthSquared, Catch::Matchers::WithinRel(25.0f));
    }

    SECTION("Dot")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        float dot = v1.dot(v2);
        REQUIRE_THAT(dot, Catch::Matchers::WithinRel(23.0f));
    }

    SECTION("Cross")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        float cross = v1.cross(v2);
        REQUIRE_THAT(cross, Catch::Matchers::WithinRel(-2.0f));
    }
}
TEST_CASE("Vec2f Min Max", "[Vec2]")
{
    SECTION("Max")
    {
        vec2<float> v(2.0f, 3.0f);
        float max = v.max();
        REQUIRE_THAT(max, Catch::Matchers::WithinRel(3.0f));
    }

    SECTION("Min")
    {
        vec2<float> v(2.0f, 3.0f);
        float min = v.min();
        REQUIRE_THAT(min, Catch::Matchers::WithinRel(2.0f));
    }

    SECTION("Max Vec2")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        vec2<float> max = v1.max(v2);
        REQUIRE_THAT(max.x, Catch::Matchers::WithinRel(4.0f));
        REQUIRE_THAT(max.y, Catch::Matchers::WithinRel(5.0f));
    }

    SECTION("Min Vec2")
    {
        vec2<float> v1(2.0f, 3.0f);
        vec2<float> v2(4.0f, 5.0f);
        vec2<float> min = v1.min(v2);
        REQUIRE_THAT(min.x, Catch::Matchers::WithinRel(2.0f));
        REQUIRE_THAT(min.y, Catch::Matchers::WithinRel(3.0f));
    }
}