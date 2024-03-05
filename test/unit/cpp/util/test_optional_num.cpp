#include "util/optional_num.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using steps::util::NoneType;
using steps::util::OptionalNum;

TEST_CASE("OptionalNum_isValue") {
    auto v3 = OptionalNum<int>(3);

    REQUIRE(v3.value() == 3);
    REQUIRE(*v3 == 3);
    REQUIRE(v3.value_or(4) == 3);
    REQUIRE(bool(v3));

    v3 = NoneType();

    REQUIRE_FALSE(v3.has_value());
    REQUIRE_FALSE(bool(v3));
    REQUIRE(v3.value_or(4) == 4);
}

TEST_CASE("OptionalNum_isNone") {
    auto n = OptionalNum<int>();
    REQUIRE_FALSE(n.has_value());
    auto n2 = OptionalNum<int>(NoneType());
    REQUIRE_FALSE(n2.has_value());
}
