#pragma once

#include <string>
#ifdef Zee_USE_TIMEMORY
#include <timemory/timemory.hpp>
#endif  // Zee_USE_TIMEMORY

namespace zee {

void timemory_cout_output(bool value);

#ifdef Zee_USE_TIMEMORY

namespace timc = tim::component;
using measurement_t = tim::component_tuple<timc::peak_rss,
                                           timc::page_rss,
                                           timc::data_rss,
                                           timc::stack_rss,
                                           timc::real_clock,
                                           timc::system_clock,
                                           timc::user_clock,
                                           timc::cpu_util>;

#else   // Zee_USE_TIMEMORY
using measurement_t = int;
#endif  // Zee_USE_TIMEMORY

class timemory_fixture {
  public:
    timemory_fixture(const std::string& name, bool store_result = false, bool auto_start = true);
    void start();
    const measurement_t& stop();

  private:
    measurement_t ti_timer;
};

}  // namespace zee
