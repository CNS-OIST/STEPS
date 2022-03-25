/*
 * Track time and memory usage for a region.
 * Print summary with min/max/avg across mpi ranks
 *
 */

#include "region_tracker.hpp"

#include <array>
#include <iomanip>
#include <utility>

#include "peak_rss.hpp"

namespace steps {
namespace util {


#ifdef STEPS_USE_MPI
void RegionTracker::init(MPI_Comm comm) {
    comm_ = comm;
}
#else
void RegionTracker::init(){}
#endif

void RegionTracker::init_region(const std::string& regionName) {
    regions_.insert({regionName, RegionTracker::tracking_t{ 0, 0, 0, 0, TimeTracker(), MemoryTracker()}});
}

void RegionTracker::start(const std::string& regionName) {
    // search for region and register if doesn't exist
    auto reg = regions_.find(regionName);
    if (reg == regions_.end()) {
        init_region(regionName);
        // retry
        start(regionName);
    }
    else {
        // start memory and time tracker (time second!)
        std::get<tracking::memTracker>(reg->second).start();
        std::get<tracking::timTracker>(reg->second).start();
    }
}

void RegionTracker::stop(const std::string& regionName) {
    // stop time and mem tracker (time first!)
    auto reg = regions_.find(regionName);
    std::get<tracking::timTracker>(reg->second).stop();
    std::get<tracking::memTracker>(reg->second).stop();

    // update count
    std::get<tracking::count>(reg->second) += 1;
    // time
    std::get<tracking::time>(reg->second) += 
        std::get<tracking::timTracker>(reg->second).diff();
    // high watermark memory (MB)
    std::get<tracking::memory>(reg->second) = 
        std::max(std::get<tracking::memory>(reg->second), peak_rss() * 1.0e-6);
    // memory delta (MB)
    std::get<tracking::memoryDelta>(reg->second) +=
        std::get<tracking::memTracker>(reg->second).diff() * 1.0e-6;

}

RegionTracker::tracking_t RegionTracker::get(const std::string& regionName) {
    auto reg = regions_.find(regionName);
    return reg->second;
}

#ifdef STEPS_USE_MPI
void RegionTracker::print(std::ostream& os) {

    int rank, n_ranks;
    MPI_Comm_rank(comm_, &rank);
    MPI_Comm_size(comm_, &n_ranks);

    for (auto reg = regions_.begin(); reg != regions_.end(); ++reg) {
        // collect stats memory
        constexpr std::size_t stat_size = 3;
        std::array<double, stat_size> data{}, data_min{}, data_max{}, data_avg{};

        // abusing the fact that data types are mapped to the lower part of the enum
        data[tracking::time] = std::get<tracking::time>(reg->second);
        data[tracking::memory] = std::get<tracking::memory>(reg->second);
        data[tracking::memoryDelta] = std::get<tracking::memoryDelta>(reg->second);

        MPI_Reduce(data.data(), data_min.data(), stat_size, MPI_DOUBLE, MPI_MIN, 0, comm_);
        MPI_Reduce(data.data(), data_max.data(), stat_size, MPI_DOUBLE, MPI_MAX, 0, comm_);
        MPI_Reduce(data.data(), data_avg.data(), stat_size, MPI_DOUBLE, MPI_SUM, 0, comm_);
        for (auto& d : data_avg) {
            d /= n_ranks;
        }

        if (rank == 0) {
            os << reg->first
                << '\t' << std::get<tracking::count>(reg->second) << " times\n"
                << "\t stats by rank                         min             max             avg\n"
                << "\t high watermark memory (MB):"
                << '\t' << std::setw(10) << data_min[tracking::memory]
                << '\t' << std::setw(10) << data_max[tracking::memory]
                << '\t' << std::setw(10) << data_avg[tracking::memory]
                << '\n'
                << "\t memory delta (MB)         :"
                << '\t' << std::setw(10) << data_min[tracking::memoryDelta]
                << '\t' << std::setw(10) << data_max[tracking::memoryDelta]
                << '\t' << std::setw(10) << data_avg[tracking::memoryDelta]
                << '\n'
                << "\t time (s)                  :"
                << '\t' << std::setw(10) << data_min[tracking::time]
                << '\t' << std::setw(10) << data_max[tracking::time]
                << '\t' << std::setw(10) << data_avg[tracking::time]
                << '\n';
        }
    }
}
#else
void RegionTracker::print(std::ostream& os) {
    for (auto reg = regions_.begin(); reg != regions_.end(); ++reg) {
        os << reg->first
            << '\t' << std::get<tracking::count>(reg->second) << " times\n"
            << "\t high watermark memory (MB):"
            << '\t' << std::setw(10) << std::get<tracking::memory>(reg->second)
            << '\n'
            << "\t memory delta (MB)         :"
            << '\t' << std::setw(10) << std::get<tracking::memoryDelta>(reg->second)
            << '\n'
            << "\t time (s)                  :"
            << '\t' << std::setw(10) << std::get<tracking::time>(reg->second)
            << '\n';
    }
}
#endif

// static members init
std::map<std::string, RegionTracker::tracking_t>
    RegionTracker::regions_{};

#ifdef STEPS_USE_MPI
MPI_Comm RegionTracker::comm_ = MPI_COMM_NULL;
#endif

} // namespace util
} // namespace steps
