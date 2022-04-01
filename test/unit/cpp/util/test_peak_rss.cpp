#include <cmath>

#include "util/tracker/peak_rss.hpp"

#include "gtest/gtest.h"

using namespace steps::util;

TEST(peak_rss, zero) {
    auto start = peak_rss();
    auto stop = peak_rss();
    ASSERT_EQ(stop - start, 0u);
}

TEST(peak_rss, MB) {
    using data_t = double;
    // page size
    // page alignement is not forced to avoid using external libraries
    constexpr std::size_t pg_s = 4096;
    // array size
    constexpr std::size_t ar_s = 3e2*pg_s;
    auto start = peak_rss();
    std::vector<data_t> v(ar_s, 1.0);
    auto stop = peak_rss();
    double measured = static_cast<double>(stop - start);
    double expected = static_cast<double>(sizeof(data_t)*ar_s);
    // not very precise, expect ~ 20 % error
    ASSERT_LE(std::abs(measured - expected)/expected, 0.2);
}

TEST(peak_rss, GB) {
    using data_t = double;
    // page size
    // page alignement is not forced to avoid using external libraries
    constexpr std::size_t pg_s = 4096;
    // array size
    constexpr std::size_t ar_s = 3e4*pg_s;
    auto start = peak_rss();
    std::vector<data_t> v(ar_s, 1.1);
    auto stop = peak_rss();
    double measured = static_cast<double>(stop - start);
    double expected = static_cast<double>(sizeof(data_t)*ar_s);
    // not very precise, expect ~ 10 % error
    ASSERT_LE(std::abs(measured - expected)/expected, 0.1);
}
