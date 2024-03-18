#include "mpi/dist/tetopsplit/kproc/event_queue.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using steps::dist::kproc::EventQueue;
using steps::dist::kproc::KProcID;

TEST_CASE("KProcEventQueue_CreateAndUpdate") {
    EventQueue queue;

    REQUIRE_THROWS_AS(queue.getEventTime(KProcID(0)), std::logic_error);
    queue.updateMaxTime(1);
    queue.update(KProcID(1), 0.5);
    REQUIRE_THAT(queue.getEventTime(KProcID(1)), Catch::Matchers::WithinULP(.5, 4));

    queue.update(KProcID(2), 0.75);
    queue.update(KProcID(3), 1.75);
    queue.update(KProcID(4), 0.2);
    queue.update(KProcID(5), 0.25);

    REQUIRE(queue.getEventTime(KProcID(4)) == 0.2);
    auto first1 = queue.getFirst();
    REQUIRE_THAT(first1.first, Catch::Matchers::WithinULP(0.2, 4));
    REQUIRE(first1.second.data() == 4);

    queue.update(KProcID(4), 0.7);
    REQUIRE_THAT(queue.getEventTime(KProcID(4)), Catch::Matchers::WithinULP(0.7, 4));
    auto first2 = queue.getFirst();
    REQUIRE_THAT(first2.first, Catch::Matchers::WithinULP(0.25, 4));
    REQUIRE(first2.second.data() == 5);

    queue.update(KProcID(1), 1.5);
    queue.update(KProcID(2), 1.75);
    queue.update(KProcID(3), 1.85);
    queue.update(KProcID(4), 1.2);
    queue.update(KProcID(5), 1.25);

    REQUIRE_THAT(queue.getEventTime(KProcID(1)), Catch::Matchers::WithinULP(1.5, 4));
    REQUIRE_THAT(queue.getEventTime(KProcID(2)), Catch::Matchers::WithinULP(1.75, 4));
    REQUIRE_THAT(queue.getEventTime(KProcID(3)), Catch::Matchers::WithinULP(1.85, 4));
    REQUIRE_THAT(queue.getEventTime(KProcID(4)), Catch::Matchers::WithinULP(1.2, 4));
    REQUIRE_THAT(queue.getEventTime(KProcID(5)), Catch::Matchers::WithinULP(1.25, 4));

    auto first3 = queue.getFirst();
    REQUIRE_THAT(first3.first,
                 Catch::Matchers::WithinULP(std::numeric_limits<double>::infinity(), 4));
    REQUIRE(first3.second.data() == 0);
}

TEST_CASE("KProcEventQueue_TimeCollision") {
    EventQueue queue;
    queue.updateMaxTime(1);
    queue.update(KProcID(2), 0.75);
    queue.update(KProcID(3), 0.75);
    queue.update(KProcID(4), 0.75);

    auto first1 = queue.getFirst();
    REQUIRE_THAT(first1.first, Catch::Matchers::WithinULP(0.75, 4));
    REQUIRE(first1.second.data() == 2);

    queue.update(KProcID(2), 0.85);

    auto first2 = queue.getFirst();
    REQUIRE_THAT(first2.first, Catch::Matchers::WithinULP(0.75, 4));
    REQUIRE(first2.second.data() == 3);

    queue.update(KProcID(3), 0.95);

    auto first3 = queue.getFirst();
    REQUIRE_THAT(first3.first, Catch::Matchers::WithinULP(0.75, 4));
    REQUIRE(first3.second.data() == 4);
}
