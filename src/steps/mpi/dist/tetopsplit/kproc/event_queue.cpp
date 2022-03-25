#include "event_queue.hpp"

#include "util/debug.hpp"
#include "util/error.hpp"

namespace steps {
namespace dist {
namespace kproc {

void EventQueue::update(const KProcID id, const osh::Real time) {
    // Remove previous event if it exists
    auto timeit = kproc_to_event_time_.find(id.data());
    if (timeit != kproc_to_event_time_.end() && timeit->second <= max_time_) {
        auto range = next_event_.equal_range(timeit->second);
        if (range.first != next_event_.end()) {
            auto it =
                std::find_if(range.first, range.second, [&id](const auto &p) {
                    return p.second == id.data();
                });
            if (it != range.second) {
                next_event_.erase(it);
            } else {
                throw std::logic_error(
                    "Can find event time but not the KProcID in the event map");
            }
        }
    }
    // Only add the event if it will happen before max_time_
    if (time != std::numeric_limits<osh::Real>::infinity() && time <= max_time_) {
        next_event_.emplace(time, id.data());
    }
    kproc_to_event_time_[id.data()] = time;
}

osh::Real EventQueue::getEventTime(const KProcID id) const {
    auto it = kproc_to_event_time_.find(id.data());
    if (it != kproc_to_event_time_.end()) {
        return it->second;
    } else {
        throw std::logic_error("Cannot find the KProcID in the event map");
    }
}

void EventQueue::updateMaxTime(const osh::Real max_time) {
    auto old_max_time = max_time_;
    if (old_max_time < max_time) {
        max_time_ = max_time;
        if(!next_event_.empty()) {
            CLOG(WARNING, "general_log") << "remain events: " << next_event_.size();
            for (auto & e : next_event_) {
                CLOG(WARNING, "general_log") << e.first << " kproc: " << e.second;
            }
            ProgErrLog("Event queue not empty!");
        }
        for(auto& kp_event : kproc_to_event_time_) {
            if(kp_event.second == std::numeric_limits<osh::Real>::infinity()) {
                continue;
            }
            if(kp_event.second <= old_max_time) {
                std::ostringstream os;
                os << "existing event time: " << kp_event.second; 
                os << " older than the previous maximum time: " << old_max_time;
                ProgErrLog(os.str());
            }
            if(kp_event.second <= max_time_) {
                next_event_.emplace(kp_event.second, kp_event.first);
            }
        }
    } else if (old_max_time > max_time) {
        ProgErrLog("Previous max time is larger than new max time, missing reset of simulation?");
    }
}

Event EventQueue::getFirst() const {
    if (next_event_.empty()) {
        return {std::numeric_limits<osh::Real>::infinity(), KProcID(0)};
    } else {
        auto event = next_event_.begin();
        // Ideally we would want to randomly select an event if several events
        // should happen at the same time. But this is rare enough that the bias
        // resulting from always selecting .begin() instead should be
        // acceptable.
        return {event->first, KProcID(event->second)};
    }
}

std::ostream &operator<<(std::ostream &ostr, const EventQueue &events) {
    return ostr << "  next_event_: " << events.next_event_
                << "\n  kproc_to_event_time_: " << events.kproc_to_event_time_;
}

} // namespace kproc
} // namespace dist
} // namespace steps
