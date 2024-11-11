#include <cmath>
#include <iostream>

#include "util/tracker/memory_tracker.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using steps::util::MemoryTracker;

TEST_CASE("MemoryTracker_zero") {
    MemoryTracker tracker{};
    tracker.start();
    tracker.stop();
    REQUIRE(tracker.diff() == 0u);
}

TEST_CASE("MemoryTracker_MB") {
    using data_t = double;
    // page size
    // page alignment is not forced to avoid using external libraries
    constexpr std::size_t pg_s = 4096;
    // array size
    constexpr std::size_t ar_s = 3e2 * pg_s;
    MemoryTracker tracker{};
    tracker.start();
    std::vector<data_t> v(ar_s, 1.0);
    std::cerr << v[0];  // prevent compiler optimization
    tracker.stop();
    auto measured = static_cast<double>(tracker.diff());
    auto expected = static_cast<double>(sizeof(data_t) * ar_s);
    // not very precise, expect ~ 20 % error
    REQUIRE(std::abs(measured - expected) / expected <= 0.27);
}

TEST_CASE("MemoryTracker_GB") {
    using data_t = double;
    // page size
    // page alignment is not forced to avoid using external libraries
    constexpr std::size_t pg_s = 4096;
    // array size
    constexpr std::size_t ar_s = 3e4 * pg_s;
    MemoryTracker tracker{};
    tracker.start();
    std::vector<data_t> v(ar_s, 1.1);
    std::cerr << v[0];  // prevent compiler optimization
    tracker.stop();
    auto measured = static_cast<double>(tracker.diff());
    auto expected = static_cast<double>(sizeof(data_t) * ar_s);
    // not very precise, expect ~ 10 % error on Linux and 25 % on Apple
#if defined(__APPLE__)
#define THRESHOLD 0.25
#else
#define THRESHOLD 0.1
#endif
    REQUIRE(std::abs(measured - expected) / expected <= THRESHOLD);
}
