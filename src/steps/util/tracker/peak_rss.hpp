/*
 * Minimal utility to get process memory usage.
 *
 */
#pragma once

#include <cstddef>

namespace steps::util {

/*
 * Extract the peak working set size (rss=resident set size).
 * ref. https://en.wikichip.org/wiki/resident_set_size
 *
 * return value is in bytes
 */

std::size_t peak_rss();


}  // namespace steps::util
