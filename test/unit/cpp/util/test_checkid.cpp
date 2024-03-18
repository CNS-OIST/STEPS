#include <string>

#include "util/checkid.hpp"
#include "util/error.hpp"

#include <catch2/catch_test_macros.hpp>

using steps::util::checkID;
using steps::util::isValidID;

std::string valid_ids[] = {"A", "z", "M", "m", "_", "a012", "B_7890", "___a", "____B", "_0Z"};

// TODO revert this temporary change for split meshes
std::string invalid_ids[] = {"0A",
                             "5z",
                             "9M",
                             "",
                             "+",
                             "@",
                             "a0!2",
                             "B_@7890",
                             "__-a",
                             "-_B",
                             /*"_.Z", */ "  ",
                             "\t",
                             "\n",
                             "a:b",
                             "#Q",
                             "3?"};

TEST_CASE("CheckID_isValidID") {
    for (const auto& s: valid_ids) {
        REQUIRE(isValidID(s));
        REQUIRE(isValidID(s.c_str()));
    }

    for (const auto& s: invalid_ids) {
        REQUIRE_FALSE(isValidID(s));
        REQUIRE_FALSE(isValidID(s.c_str()));
    }
}

TEST_CASE("CheckID_checkID") {
    for (const auto& s: valid_ids) {
        REQUIRE_NOTHROW(checkID(s));
        REQUIRE_NOTHROW(checkID(s.c_str()));
    }

    for (const auto& s: invalid_ids) {
        REQUIRE_THROWS_AS(checkID(s), steps::ArgErr);
        REQUIRE_THROWS_AS(checkID(s.c_str()), steps::ArgErr);
    }
}
