/*
 * Track maximum memory usage between 2 points
 *
 */

#include "memory_tracker.hpp"
#include "peak_rss.hpp"

namespace steps::util {

void MemoryTracker::start() {
    init_ = peak_rss();
}

void MemoryTracker::stop() {
    final_ = peak_rss();
}

/*
 * return value is in bytes
 */
std::size_t MemoryTracker::diff() const {
    return final_ - init_;
}

}  // namespace steps::util
