#include "steps/util/strong_id.hpp"

#include <cstdint>

#include "gtest/gtest.h"

using namespace steps::util;

class strong_idTest : public ::testing::Test {
public:
    static constexpr struct random_id_trait_t {} random_id_trait {};
    static constexpr struct not_a_number_t {} not_a_number {};
    static constexpr struct also_not_a_number_t {} also_not_a_number {};
};

TEST_F(strong_idTest, is_promotion_SameFamilyType) {
    ASSERT_TRUE((is_promotion<char, int>::value));
    ASSERT_TRUE((is_promotion<long long, std::int64_t>::value));
    ASSERT_TRUE((is_promotion<float, double>::value));

    ASSERT_FALSE((is_promotion<double, float>::value));
    ASSERT_FALSE((is_promotion<std::uint64_t, std::int64_t>::value));
    ASSERT_FALSE((is_promotion<int, char>::value));

    // not sure if this should fail or not, if it should pass we need to add
    // is_integral, is_arithmetic, and std::numeric_limits<>::is_specialized
    // in addition to numeric_limits<>::max() (GC)
    ASSERT_FALSE((is_promotion<strong_id<long, random_id_trait_t>, std::int64_t>::value));
}

TEST_F(strong_idTest, is_promotion_MixedFamilyType) {
    ASSERT_FALSE((is_promotion<float, long>::value));
    ASSERT_FALSE((is_promotion<int, double>::value));
    ASSERT_FALSE((is_promotion<not_a_number_t, double>::value));
    ASSERT_FALSE((is_promotion<not_a_number_t, also_not_a_number_t>::value));
}

TEST_F(strong_idTest, numeric_limits) {
    using atype = long;
    ASSERT_EQ((std::numeric_limits<strong_id<atype, random_id_trait_t>>::max()),
        std::numeric_limits<atype>::max());
}
