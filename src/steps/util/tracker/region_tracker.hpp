/*
 * Track time and memory usage for a region.
 * Print summary with min/max/avg across mpi ranks
 *
 */
#pragma once

#include "memory_tracker.hpp"
#include "time_tracker.hpp"

#include <cstddef>
#include <cstdint>
#include <iostream>
#include <map>
#ifdef STEPS_USE_MPI
#include <mpi.h>
#endif
#include <tuple>

namespace steps {
namespace util {


class RegionTracker {

    using tracking_t = std::tuple<double, double, double, std::uint32_t, TimeTracker, MemoryTracker>;
private:
    enum tracking {time, memory, memoryDelta, count, timTracker, memTracker};

    // region name -> tracking tuple
    static std::map<std::string, tracking_t> regions_;

#ifdef STEPS_USE_MPI
    // MPI communicator
    static MPI_Comm comm_;
#endif

public:

#ifdef STEPS_USE_MPI
    static void init(MPI_Comm comm = MPI_COMM_WORLD);
#else
    static void init();
#endif
    static void init_region(const std::string& name);
    static void start(const std::string& name);
    static void stop(const std::string& name);
    static tracking_t get(const std::string& name);
    static void print(std::ostream& os = std::cout);

};

} // namespace util
} // namespace steps
