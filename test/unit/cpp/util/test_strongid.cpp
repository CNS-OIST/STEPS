#include "steps/util/strong_id.hpp"
#include "steps/util/strong_ra.hpp"

#include <numeric>

#include "gtest/gtest.h"

using steps::util::is_promotion;
using steps::util::strong_id;

class strong_idTest: public ::testing::Test {
  public:
    static constexpr struct random_id_trait_t {
    } random_id_trait{};
    static constexpr struct not_a_number_t {
    } not_a_number{};
    static constexpr struct also_not_a_number_t { } also_not_a_number{}; };


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

TEST_F(strong_idTest, size) {
    static_assert(sizeof(uint64_t) == sizeof(strong_id<uint64_t, random_id_trait_t>),
                  "struct size should be the same than embedded integral type");
    static_assert(sizeof(int32_t) == sizeof(strong_id<int32_t, random_id_trait_t>),
                  "struct size should be the same than embedded integral type");
}

TEST_F(strong_idTest, constructor) {
    {
        // default constructor
        strong_id<long, random_id_trait_t> id;
        ASSERT_TRUE(id.get() == decltype(id)::unknown_value());
        ASSERT_TRUE(id.unknown());
        ASSERT_FALSE(id.valid());
    }
    {
        // uniform initialization
        strong_id<long, random_id_trait_t> id{};
        ASSERT_TRUE(id.get() == decltype(id)::unknown_value());
        ASSERT_TRUE(id.unknown());
        ASSERT_FALSE(id.valid());
    }

    {
        // assign valid values 0 and 42
        using sid = strong_id<long, random_id_trait_t>;
        sid id;
        for (auto val: {sid(0), sid(42)}) {
            id = val;
            ASSERT_EQ(id.get(), val);
            ASSERT_FALSE(id.unknown());
            ASSERT_TRUE(id.valid());
        }
    }
}

TEST_F(strong_idTest, vector) {
    using sid = strong_id<long, random_id_trait_t>;
    using vec = std::vector<sid>;
    vec v(1);
    ASSERT_EQ(v.back(), sid::unknown_value());
    v.resize(2);
    ASSERT_EQ(v.back(), sid::unknown_value());
    v.emplace_back();
    ASSERT_EQ(v.back(), sid::unknown_value());
}

TEST_F(strong_idTest, array) {
    using sid = strong_id<long, random_id_trait_t>;
    using array = std::array<sid, 2>;
    {
        array a{std::nullopt, std::nullopt};
        ASSERT_EQ(a[0].get(), sid::unknown_value());
        ASSERT_EQ(a[1].get(), sid::unknown_value());
    }
    {
        array a{{{}, {}}};
        ASSERT_EQ(a[0].get(), sid::unknown_value());
        ASSERT_EQ(a[1].get(), sid::unknown_value());
    }
}

TEST(strong_ra, vector_ref) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    std::vector<int> data({1, 2, 3});

    auto sdata = steps::util::make_strong_random_accessor<sid>(data);
    static_assert(!sdata.container_owned, "vector should not be copied");
    ASSERT_EQ(&data, &sdata.container());

    ASSERT_EQ(data.size(), sdata.size());
    ASSERT_EQ(sdata[sid(2)], 3);
    data[2] = 4;
    ASSERT_EQ(sdata[sid(2)], 4);
}

TEST(strong_ra, vector_const_ref) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    const std::vector<int> data({1, 2, 3});

    auto sdata = steps::util::make_strong_random_accessor<sid>(data);
    static_assert(!sdata.container_owned, "vector should not be copied");
    ASSERT_EQ(&data, &sdata.container());
    ASSERT_EQ(data.size(), sdata.size());
    ASSERT_EQ(sdata[sid(2)], 3);
    ASSERT_EQ(sdata.at(sid(2)), 3);
    ASSERT_EQ(sdata.front(), 1);
    ASSERT_EQ(sdata.back(), 3);
    ASSERT_NE(sdata.data(), nullptr);
    ASSERT_EQ(std::accumulate(sdata.begin(), sdata.end(), 0), 6);
}

TEST(strong_ra, vector_rvalue_ref) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    std::vector<int> data({1, 2, 3});
    auto data_ptr = data.data();
    auto sdata = steps::util::make_strong_random_accessor<sid>(std::move(data));
    static_assert(sdata.container_owned, "vector should be owned by the strong_ra instance");

    ASSERT_EQ(sdata.size(), 3);
    // vector must have been moved
    ASSERT_EQ(sdata.data(), data_ptr);
    ASSERT_EQ(data.size(), 0);
    // vector should be modifiable
    sdata.front() = 0;
}

TEST(strong_ra, vector_instantiation) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    {
        auto strong_vector = steps::util::make_strong_random_accessor<sid, std::vector<int>>(3);
        static_assert(strong_vector.container_owned, "vector should be owned");
        ASSERT_EQ(strong_vector.size(), 3);
    }
    {
        auto strong_vector = steps::util::make_strong_random_accessor<sid, std::vector<int>>(3, 5);
        static_assert(strong_vector.container_owned, "vector should be owned");
        ASSERT_EQ(strong_vector.size(), 3);
        ASSERT_EQ(strong_vector.back(), 5);
    }
}

TEST(strong_ra, span) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    std::vector<int> data({1, 2, 3});
    gsl::span<int> span_data{data.data(), data.size()};
    auto strong_data = steps::util::make_strong_random_accessor<sid>(span_data);
    ASSERT_EQ(data.size(), strong_data.size());
    ASSERT_EQ(data.data(), strong_data.data());
    static_assert(strong_data.container_owned, "vector should not be copied");
}

TEST(strong_ra, create_no_arg) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;

    auto strong_vector = steps::util::strong_random_access<sid, std::vector<double>>();
    static_assert(strong_vector.container_owned, "vector should be owned");
    strong_vector.container().push_back(42.4);
    ASSERT_EQ(strong_vector.size(), 1);
}

#ifdef STEPS_USE_DIST_MESH

TEST(strong_ra, omega_h_write) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    Omega_h::Write<Omega_h::LO> write(10, 3);
    auto strong_vector = steps::util::make_strong_random_accessor<sid>(write);
    ASSERT_EQ(strong_vector.size(), 10);
    ASSERT_EQ(strong_vector[sid(0)], 3);
    strong_vector[sid(0)] = 0;
    ASSERT_EQ(strong_vector[sid(0)], 0);
    ASSERT_EQ(write[0], 0);
}

TEST(strong_ra, omega_h_read) {
    using sid = strong_id<long, strong_idTest::random_id_trait_t>;
    Omega_h::Write<Omega_h::LO> write(10, 3);
    Omega_h::Read<Omega_h::LO> read(write);
    auto strong_vector = steps::util::make_strong_random_accessor<sid>(read);
    ASSERT_EQ(strong_vector.size(), 10);
    ASSERT_EQ(strong_vector[sid(0)], 3);
}

#endif  // STEPS_USE_DIST_MESH
