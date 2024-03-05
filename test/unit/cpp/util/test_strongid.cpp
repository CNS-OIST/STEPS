#include "steps/util/strong_id.hpp"
#include "steps/util/strong_ra.hpp"

#include <numeric>

#include <catch2/catch_test_macros.hpp>

using steps::util::is_promotion;
using steps::util::strong_id;

static constexpr struct random_id_trait_t {
} random_id_trait{};
static constexpr struct not_a_number_t {
} not_a_number{};
static constexpr struct also_not_a_number_t {
} also_not_a_number{};


TEST_CASE("strong_idTest_is_promotion_SameFamilyType") {
    REQUIRE((is_promotion<char, int>::value));
    REQUIRE((is_promotion<long long, std::int64_t>::value));
    REQUIRE((is_promotion<float, double>::value));

    REQUIRE_FALSE((is_promotion<double, float>::value));
    REQUIRE_FALSE((is_promotion<std::uint64_t, std::int64_t>::value));
    REQUIRE_FALSE((is_promotion<int, char>::value));

    // not sure if this should fail or not, if it should pass we need to add
    // is_integral, is_arithmetic, and std::numeric_limits<>::is_specialized
    // in addition to numeric_limits<>::max() (GC)
    REQUIRE_FALSE((is_promotion<strong_id<long, random_id_trait_t>, std::int64_t>::value));
}

TEST_CASE("strong_idTest_is_promotion_MixedFamilyType") {
    REQUIRE_FALSE((is_promotion<float, long>::value));
    REQUIRE_FALSE((is_promotion<int, double>::value));
    REQUIRE_FALSE((is_promotion<not_a_number_t, double>::value));
    REQUIRE_FALSE((is_promotion<not_a_number_t, also_not_a_number_t>::value));
}

TEST_CASE("strong_idTest_numeric_limits") {
    using atype = long;
    REQUIRE((std::numeric_limits<strong_id<atype, random_id_trait_t>>::max()) ==
            std::numeric_limits<atype>::max());
}

TEST_CASE("strong_idTest_size") {
    static_assert(sizeof(uint64_t) == sizeof(strong_id<uint64_t, random_id_trait_t>),
                  "struct size should be the same than embedded integral type");
    static_assert(sizeof(int32_t) == sizeof(strong_id<int32_t, random_id_trait_t>),
                  "struct size should be the same than embedded integral type");
}

TEST_CASE("strong_idTest_constructor") {
    {
        // default constructor
        strong_id<long, random_id_trait_t> id;
        REQUIRE(id.get() == decltype(id)::unknown_value());
        REQUIRE(id.unknown());
        REQUIRE_FALSE(id.valid());
    }
    {
        // uniform initialization
        strong_id<long, random_id_trait_t> id{};
        REQUIRE(id.get() == decltype(id)::unknown_value());
        REQUIRE(id.unknown());
        REQUIRE_FALSE(id.valid());
    }

    {
        // assign valid values 0 and 42
        using sid = strong_id<long, random_id_trait_t>;
        sid id;
        for (auto val: {sid(0), sid(42)}) {
            id = val;
            REQUIRE(id.get() == val);
            REQUIRE_FALSE(id.unknown());
            REQUIRE(id.valid());
        }
    }
}

TEST_CASE("strong_idTest_vector") {
    using sid = strong_id<long, random_id_trait_t>;
    using vec = std::vector<sid>;
    vec v(1);
    REQUIRE(v.back() == sid::unknown_value());
    v.resize(2);
    REQUIRE(v.back() == sid::unknown_value());
    v.emplace_back();
    REQUIRE(v.back() == sid::unknown_value());
}

TEST_CASE("strong_idTest_array") {
    using sid = strong_id<long, random_id_trait_t>;
    using array = std::array<sid, 2>;
    {
        array a{std::nullopt, std::nullopt};
        REQUIRE(a[0].get() == sid::unknown_value());
        REQUIRE(a[1].get() == sid::unknown_value());
    }
    {
        array a{{{}, {}}};
        REQUIRE(a[0].get() == sid::unknown_value());
        REQUIRE(a[1].get() == sid::unknown_value());
    }
}

TEST_CASE("strong_ra_vector_ref") {
    using sid = strong_id<long, random_id_trait_t>;
    std::vector<int> data({1, 2, 3});

    auto sdata = steps::util::make_strong_random_accessor<sid>(data);
    static_assert(!sdata.container_owned, "vector should not be copied");
    REQUIRE(&data == &sdata.container());

    REQUIRE(data.size() == sdata.size());
    REQUIRE(sdata[sid(2)] == 3);
    data[2] = 4;
    REQUIRE(sdata[sid(2)] == 4);
}

TEST_CASE("strong_ra_vector_const_ref") {
    using sid = strong_id<long, random_id_trait_t>;
    const std::vector<int> data({1, 2, 3});

    auto sdata = steps::util::make_strong_random_accessor<sid>(data);
    static_assert(!sdata.container_owned, "vector should not be copied");
    REQUIRE(&data == &sdata.container());
    REQUIRE(data.size() == sdata.size());
    REQUIRE(sdata[sid(2)] == 3);
    REQUIRE(sdata.at(sid(2)) == 3);
    REQUIRE(sdata.front() == 1);
    REQUIRE(sdata.back() == 3);
    REQUIRE(sdata.data() != nullptr);
    REQUIRE(std::accumulate(sdata.begin(), sdata.end(), 0) == 6);
}

TEST_CASE("strong_ra_vector_rvalue_ref") {
    using sid = strong_id<long, random_id_trait_t>;
    std::vector<int> data({1, 2, 3});
    auto data_ptr = data.data();
    auto sdata = steps::util::make_strong_random_accessor<sid>(std::move(data));
    static_assert(sdata.container_owned, "vector should be owned by the strong_ra instance");

    REQUIRE(sdata.size() == 3);
    // vector must have been moved
    REQUIRE(sdata.data() == data_ptr);
    REQUIRE(data.size() == 0);
    // vector should be modifiable
    sdata.front() = 0;
}

TEST_CASE("strong_ra_vector_instantiation") {
    using sid = strong_id<long, random_id_trait_t>;
    {
        auto strong_vector = steps::util::make_strong_random_accessor<sid, std::vector<int>>(3);
        static_assert(strong_vector.container_owned, "vector should be owned");
        REQUIRE(strong_vector.size() == 3);
    }
    {
        auto strong_vector = steps::util::make_strong_random_accessor<sid, std::vector<int>>(3, 5);
        static_assert(strong_vector.container_owned, "vector should be owned");
        REQUIRE(strong_vector.size() == 3);
        REQUIRE(strong_vector.back() == 5);
    }
}

TEST_CASE("strong_ra_span") {
    using sid = strong_id<long, random_id_trait_t>;
    std::vector<int> data({1, 2, 3});
    gsl::span<int> span_data{data.data(), data.size()};
    auto strong_data = steps::util::make_strong_random_accessor<sid>(span_data);
    REQUIRE(data.size() == strong_data.size());
    REQUIRE(data.data() == strong_data.data());
    static_assert(strong_data.container_owned, "vector should not be copied");
}

TEST_CASE("strong_ra_create_no_arg") {
    using sid = strong_id<long, random_id_trait_t>;

    auto strong_vector = steps::util::strong_random_access<sid, std::vector<double>>();
    static_assert(strong_vector.container_owned, "vector should be owned");
    strong_vector.container().push_back(42.4);
    REQUIRE(strong_vector.size() == 1);
}

#ifdef STEPS_USE_DIST_MESH

TEST_CASE("strong_ra_omega_h_write") {
    using sid = strong_id<long, random_id_trait_t>;
    Omega_h::Write<Omega_h::LO> write(10, 3);
    auto strong_vector = steps::util::make_strong_random_accessor<sid>(write);
    REQUIRE(strong_vector.size() == 10);
    REQUIRE(strong_vector[sid(0)] == 3);
    strong_vector[sid(0)] = 0;
    REQUIRE(strong_vector[sid(0)] == 0);
    REQUIRE(write[0] == 0);
}

TEST_CASE("strong_ra_omega_h_read") {
    using sid = strong_id<long, random_id_trait_t>;
    Omega_h::Write<Omega_h::LO> write(10, 3);
    Omega_h::Read<Omega_h::LO> read(write);
    auto strong_vector = steps::util::make_strong_random_accessor<sid>(read);
    REQUIRE(strong_vector.size() == 10);
    REQUIRE(strong_vector[sid(0)] == 3);
}

#endif  // STEPS_USE_DIST_MESH
