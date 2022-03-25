/*
 * Track time usage between 2 points
 *
 */
#pragma once

#include <chrono>

namespace steps {
namespace util {

/*
 * return value is in seconds
 */

class TimeTracker {

public:
    void start();
    void stop();
    double diff();
private:
    std::chrono::steady_clock::time_point init_{};
    std::chrono::steady_clock::time_point final_{};
};

} // namespace util
} // namespace steps
