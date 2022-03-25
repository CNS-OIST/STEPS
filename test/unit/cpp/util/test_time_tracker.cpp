#include <cmath>
#include <thread>
#include <chrono>

#include "util/tracker/time_tracker.hpp"

#include "gtest/gtest.h"

using namespace steps::util;

TEST(TimeTracker, zero) {
    TimeTracker tracker{};
    tracker.start();
    tracker.stop();
    ASSERT_EQ(tracker.diff(), 0);
}

// expected good resolution in the order seconds

TEST(TimeTracker, s) {
    using namespace std::chrono_literals;
    TimeTracker tracker{};
    tracker.start();
    std::this_thread::sleep_for(1000ms);
    tracker.stop();
    double measured = tracker.diff();
    double expected = 1.0;
    ASSERT_LE(std::abs(measured - expected)/expected, 0.05);
}
