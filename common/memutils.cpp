#include "memutils.hpp"

namespace zee {

#ifdef Zee_USE_TIMEMORY
void timemory_cout_output(bool value) {
    tim::settings::cout_output() = value;
}
#else
void timemory_cout_output(bool) {}
#endif  // ZEE_USE_TIMEMORY

timemory_fixture::timemory_fixture(const std::string& name, bool store_result, bool auto_start)
#ifdef Zee_USE_TIMEMORY
    : ti_timer(name, store_result) {
    if (auto_start) {
        start();
    }
}
#else
    : ti_timer(0) {
    static_cast<void>(name);
    static_cast<void>(store_result);
    static_cast<void>(auto_start);
}
#endif  // Zee_USE_TIMEMORY

void timemory_fixture::start() {
#ifdef Zee_USE_TIMEMORY
    ti_timer.start();
#endif  // Zee_USE_TIMEMORY
}


const measurement_t& timemory_fixture::stop() {
#ifdef Zee_USE_TIMEMORY
    ti_timer.stop();
#endif
    return ti_timer;
}

}  // namespace zee
