#include <string>

#include "steps/error.hpp"
#include "steps/util/checkid.hpp"

#include "gtest/gtest.h"

using namespace steps::util;

std::string valid_ids[]={
    "A", "z", "M", "m", "_",
    "a012", "B_7890", "___a", "____B", "_0Z"
};

std::string invalid_ids[]={
    "0A", "5z", "9M", "", "+", "@", 
    "a0!2", "B_@7890", "__-a", "-_B", "_.Z", "  ",
    "\t", "\n", "a:b", "#Q", "3?"
};

TEST(CheckID,isValidID) {
    for (const auto &s: valid_ids) {
        ASSERT_TRUE(isValidID(s));
        ASSERT_TRUE(isValidID(s.c_str()));
    }

    for (const auto &s: invalid_ids) {
        ASSERT_FALSE(isValidID(s));
        ASSERT_FALSE(isValidID(s.c_str()));
    }
}

TEST(CheckID,checkID) {
    for (const auto &s: valid_ids) {
        EXPECT_NO_THROW(checkID(s));
        EXPECT_NO_THROW(checkID(s.c_str()));
    }

    for (const auto &s: invalid_ids) {
        ASSERT_THROW(checkID(s),steps::ArgErr);
        ASSERT_THROW(checkID(s.c_str()),steps::ArgErr);
    }
}
