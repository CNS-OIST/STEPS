#include "util/flat_multimap.hpp"

#include <catch2/catch_approx.hpp>
#include <catch2/catch_template_test_macros.hpp>

using namespace Catch::literals;  // NOLINT
#include <catch2/matchers/catch_matchers_floating_point.hpp>


using steps::fmm_ab2c_padding;
using steps::fmm_stl;

template <int Mask>
struct Policy {
    static constexpr int value = Mask;
};

#ifdef STEPS_USE_DIST_MESH
using steps::fmm_osh;
#define FMM_BACKENDS steps::util::OSH, steps::util::STL
typedef std::tuple<Policy<fmm_osh>,
                   Policy<fmm_osh | fmm_ab2c_padding>,
                   Policy<fmm_stl>,
                   Policy<fmm_stl | fmm_ab2c_padding>>
    fmm_test_policies;
#else
#define FMM_BACKENDS steps::util::STL
using fmm_test_policies = std::tuple<Policy<fmm_stl>, Policy<fmm_stl | fmm_ab2c_padding>>;
#endif  // STEPS_USE_DIST_MESH

TEMPLATE_LIST_TEST_CASE("size 1", "[vds]", fmm_test_policies) {
    steps::util::flat_multimap<double, 1, TestType::value> v({1, 3, 0, 2}, 3.);

    REQUIRE(v.size() == 4);
    REQUIRE(v.size(0) == 1);
    REQUIRE(v.size(1) == 3);
    REQUIRE(v.size(2) == 0);
    REQUIRE(v.size(3) == 2);
    REQUIRE(v.empty() == false);

    REQUIRE(v(0, 0) == 3.);
    REQUIRE(v(1, 0) == 3.);
    REQUIRE(v(1, 1) == 3.);
    REQUIRE(v(1, 2) == 3.);
    REQUIRE(v(3, 0) == 3.);
    REQUIRE(v(3, 1) == 3.);

    REQUIRE(v.empty(0) == false);
    REQUIRE(v.empty(1) == false);
    REQUIRE(v.empty(2) == true);
    REQUIRE(v.empty(3) == false);

    REQUIRE(v.at(0).size() == 1);
    REQUIRE(v.at(1).size() == 3);
    REQUIRE(v.at(2).size() == 0);
    REQUIRE(v.at(3).size() == 2);

    {
        double value = 1.;
        v(0, 0) = value++;
        v(1, 0) = value++;
        v(1, 1) = value++;
        v(1, 2) = value++;
        v(3, 0) = value++;
        v(3, 1) = value++;
    }

    REQUIRE(v(0, 0) == 1.);
    REQUIRE(v(1, 0) == 2.);
    REQUIRE(v(1, 1) == 3.);
    REQUIRE(v(1, 2) == 4.);
    REQUIRE(v(3, 0) == 5.);
    REQUIRE(v(3, 1) == 6.);

    {
        // use constant iterators
        size_t num_items{};
        double sum{};
        for (const auto& range: v) {
            for (const auto value: range) {
                ++num_items;
                sum += value;
            }
        }
        REQUIRE(num_items == 6);
        REQUIRE(static_cast<int>(sum) == 21);
    }
    {
        // use iterators to increment all values
        for (auto& range: v) {
            for (auto& value: range) {
                ++value;
            }
        }
        size_t num_items{};
        double sum{};
        for (const auto& range: v) {
            for (auto value: range) {
                sum += value;
                ++num_items;
            }
        }
        REQUIRE(num_items == 6);
        REQUIRE(static_cast<int>(sum) == 27);
    }

    {
        v.assign(-3.);
        REQUIRE(v(1, 0) == -3._a);
        REQUIRE(v(1, 1) == -3._a);
        REQUIRE(v(1, 2) == -3._a);
        REQUIRE(v(3, 0) == -3._a);
    }
}

TEMPLATE_LIST_TEST_CASE("size 3", "[vds]", fmm_test_policies) {
    using fmm_type = typename steps::util::flat_multimap<double, 3, TestType::value>;
    fmm_type v({1, 3, 0, 2}, 3.);

    REQUIRE(v.size() == 4);
    REQUIRE(v.size(0) == 1);
    REQUIRE(v.size(1) == 3);
    REQUIRE(v.size(2) == 0);
    REQUIRE(v.size(3) == 2);
    REQUIRE(v.empty() == false);

    REQUIRE(v(0, 0).size() == 3);
    REQUIRE(v(0, 0)[0] == 3.);
    REQUIRE(v(0, 0)[1] == 3.);
    REQUIRE(v(0, 0)[2] == 3.);

    REQUIRE(v(1, 0).size() == 3);
    REQUIRE(v(1, 0)[0] == 3.);
    REQUIRE(v(1, 0)[1] == 3.);
    REQUIRE(v(1, 0)[2] == 3.);

    REQUIRE(v(1, 1).size() == 3);
    REQUIRE(v(1, 1)[0] == 3.);
    REQUIRE(v(1, 1)[1] == 3.);
    REQUIRE(v(1, 1)[2] == 3.);

    REQUIRE(v(1, 2).size() == 3);
    REQUIRE(v(1, 2)[0] == 3.);
    REQUIRE(v(1, 2)[1] == 3.);
    REQUIRE(v(1, 2)[2] == 3.);

    REQUIRE(v(3, 0).size() == 3);
    REQUIRE(v(3, 0)[0] == 3.);
    REQUIRE(v(3, 0)[1] == 3.);
    REQUIRE(v(3, 0)[2] == 3.);

    REQUIRE(v(3, 1).size() == 3);
    REQUIRE(v(3, 1)[0] == 3.);
    REQUIRE(v(3, 1)[1] == 3.);
    REQUIRE(v(3, 1)[2] == 3.);

    {
        double value = 7;
        size_t num_elements{};
        size_t num_data{};
        size_t num_items{};
        for (auto& element: v) {
            ++num_elements;
            for (auto& data: element) {
                ++num_data;
                for (auto& item: data) {
                    item = value++;
                    ++num_items;
                }
            }
        }
        REQUIRE(num_elements == 4);
        REQUIRE(num_data == 6);
        REQUIRE(num_items == 18);
    }
    {
        REQUIRE(v(0, 0).size() == 3);
        REQUIRE(v(0, 0)[0] == 7.);
        REQUIRE(v(0, 0)[1] == 8.);
        REQUIRE(v(0, 0)[2] == 9.);

        REQUIRE(v(1, 0).size() == 3);
        REQUIRE(v(1, 0)[0] == 10.);
        REQUIRE(v(1, 0)[1] == 11.);
        REQUIRE(v(1, 0)[2] == 12.);

        REQUIRE(v(1, 1).size() == 3);
        REQUIRE(v(1, 1)[0] == 13.);
        REQUIRE(v(1, 1)[1] == 14.);
        REQUIRE(v(1, 1)[2] == 15.);

        REQUIRE(v(1, 2).size() == 3);
        REQUIRE(v(1, 2)[0] == 16.);
        REQUIRE(v(1, 2)[1] == 17.);
        REQUIRE(v(1, 2)[2] == 18.);

        REQUIRE(v(3, 0).size() == 3);
        REQUIRE(v(3, 0)[0] == 19.);
        REQUIRE(v(3, 0)[1] == 20.);
        REQUIRE(v(3, 0)[2] == 21.);

        REQUIRE(v(3, 1).size() == 3);
        REQUIRE(v(3, 1)[0] == 22.);
        REQUIRE(v(3, 1)[1] == 23.);
        REQUIRE(v(3, 1)[2] == 24.);
    }
    {
        auto it = v.begin();
        auto begin = it++;
        std::advance(it, 2);
        REQUIRE(begin == v.begin());
        REQUIRE(++it == v.end());
        REQUIRE(it == v.end());
        REQUIRE(it++ == v.end());
    }
    {
        // test copy constructor;
        fmm_type v2(v);
        REQUIRE(v2.size() == 4);
        REQUIRE(v2(0, 0).size() == 3);
    }
    {
        // test assignment operator
        fmm_type v2;
        v2 = std::move(v);
        REQUIRE(v2.size() == 4);
        REQUIRE(v2(0, 0).size() == 3);
    }
}
