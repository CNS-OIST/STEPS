#include <chrono>
#include <cmath>
#include <thread>

#include "util/tracker/time_tracker.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

TEST_CASE("TimeTracker_zero") {
    steps::util::TimeTracker tracker{};
    tracker.start();
    tracker.stop();
    REQUIRE_THAT(tracker.diff(), Catch::Matchers::WithinULP(0., 4));
}

// expected good resolution in the order seconds

TEST_CASE("TimeTracker_s") {
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
    REQUIRE(std::abs(measured - expected) / expected <= THRESHOLD);
}
