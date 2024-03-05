/*
 * Track time usage between 2 points.
 *
 */

#include "time_tracker.hpp"

namespace steps::util {

void TimeTracker::start() {
    init_ = std::chrono::steady_clock::now();
}

void TimeTracker::stop() {
    final_ = std::chrono::steady_clock::now();
}

/*
 * return value is in seconds
 */

double TimeTracker::diff() {
    return static_cast<double>(
               std::chrono::duration_cast<std::chrono::microseconds>(final_ - init_).count()) *
           1e-6;
}

}  // namespace steps::util
