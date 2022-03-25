#include "mol_state.hpp"
#include "util/debug.hpp"

namespace steps {
namespace dist {

std::ostream& operator<<(std::ostream& os, const Occupancy& o) {
    os << "correction_: " << o.corrections_ << '\n';
    os << "ids_: " << o.ids_ << '\n';
    os << "start_time_: " << o.start_time_ << '\n';
    return os;
}

}  // namespace dist
}  // namespace steps
