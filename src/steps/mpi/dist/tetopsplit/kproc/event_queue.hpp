#pragma once

#include <iosfwd>
#include <unordered_map>

#include "kproc_id.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"

namespace steps {
namespace dist {
namespace kproc {

using EventTime = osh::Real;
using Event = std::pair<EventTime, kproc::KProcID>;

/**
 * \brief A small utility struct to keeep track of the times at which
 * future events should happen.
 */
class EventQueue {
  public:
    EventQueue() = default;

    void clear() {
        next_event_.clear();
        kproc_to_event_time_.clear();
        max_time_ = -std::numeric_limits<osh::Real>::epsilon();
    }

    void update(const KProcID id, const osh::Real time);

    osh::Real getEventTime(const KProcID id) const;
    void updateMaxTime(const osh::Real max_time);
    Event getFirst() const;

    friend std::ostream& operator<<(std::ostream& ostr, const EventQueue& events);

  private:
    std::multimap<osh::Real, unsigned> next_event_ {};
    std::unordered_map<unsigned, osh::Real> kproc_to_event_time_ {};

    // maximum time that an event is put into the multimap
    osh::Real max_time_ {-std::numeric_limits<osh::Real>::epsilon()};
};

/**
 * \a EventQueue pretty printer
 */
std::ostream &operator<<(std::ostream &ostr, const EventQueue &events);

} // namespace kproc
} // namespace dist
} // namespace steps
