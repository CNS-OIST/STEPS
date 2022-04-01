/*
 * Track memory usage between 2 points
 *
 */
#pragma once

#include <cstddef>

namespace steps {
namespace util {


class MemoryTracker {

public:
    void start();
    void stop();

    /*
     * return value is in bytes
     */
    std::size_t diff();
private:
    std::size_t init_{};
    std::size_t final_{};
};

} // namespace util
} // namespace steps
