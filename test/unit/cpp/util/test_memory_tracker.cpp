#include <cmath>

#include "util/tracker/memory_tracker.hpp"

#include "gtest/gtest.h"

using steps::util::MemoryTracker;

TEST(MemoryTracker, zero) {
    MemoryTracker tracker{};
    tracker.start();
    tracker.stop();
    ASSERT_EQ(tracker.diff(), 0u);
}

TEST(MemoryTracker, MB) {
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
    double measured = static_cast<double>(tracker.diff());
    double expected = static_cast<double>(sizeof(data_t) * ar_s);
    // not very precise, expect ~ 20 % error
    ASSERT_LE(std::abs(measured - expected) / expected, 0.2);
}

TEST(MemoryTracker, GB) {
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
    double measured = static_cast<double>(tracker.diff());
    double expected = static_cast<double>(sizeof(data_t) * ar_s);
    // not very precise, expect ~ 10 % error on Linux and 25 % on Apple
#if defined(__APPLE__)
#define THRESHOLD 0.25
#else
#define THRESHOLD 0.1
#endif
    ASSERT_LE(std::abs(measured - expected) / expected, THRESHOLD);
}
