/*
 * Track memory usage between 2 points
 *
 */
#pragma once

#include <cstddef>

namespace steps::util {


class MemoryTracker {
  public:
    void start();
    void stop();

    /*
     * return value is in bytes
     */
    std::size_t diff() const;

  private:
    std::size_t init_{};
    std::size_t final_{};
};

}  // namespace steps::util
