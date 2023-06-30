#include <chrono>
#include <cmath>
#include <thread>

#include "util/tracker/time_tracker.hpp"

#include "gtest/gtest.h"

TEST(TimeTracker, zero) {
    steps::util::TimeTracker tracker{};
    tracker.start();
    tracker.stop();
    ASSERT_EQ(tracker.diff(), 0);
}

// expected good resolution in the order seconds

TEST(TimeTracker, s) {
    using namespace std::chrono_literals;  // NOLINT
    steps::util::TimeTracker tracker{};
    tracker.start();
    std::this_thread::sleep_for(1000ms);
    tracker.stop();
    double measured = tracker.diff();
    double expected = 1.0;
#if defined(__APPLE__)
#define THRESHOLD 0.2
#else
#define THRESHOLD 0.05
#endif
    ASSERT_LE(std::abs(measured - expected) / expected, THRESHOLD);
}
