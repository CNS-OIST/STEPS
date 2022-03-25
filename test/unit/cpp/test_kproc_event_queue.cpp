#include "mpi/dist/tetopsplit/kproc/event_queue.hpp"

#include "gtest/gtest.h"

using namespace steps::dist::kproc;

TEST(KProcEventQueue, CreateAndUpdate) {
    EventQueue queue;

    ASSERT_THROW(queue.getEventTime(KProcID(0)), std::logic_error);
    queue.updateMaxTime(1);
    queue.update(KProcID(1), 0.5);
    ASSERT_EQ(queue.getEventTime(KProcID(1)), .5);

    queue.update(KProcID(2), 0.75);
    queue.update(KProcID(3), 1.75);
    queue.update(KProcID(4), 0.2);
    queue.update(KProcID(5), 0.25);

    ASSERT_EQ(queue.getEventTime(KProcID(4)), 0.2);
    auto first1 = queue.getFirst();
    ASSERT_EQ(first1.first, 0.2);
    ASSERT_EQ(first1.second.data(), 4);

    queue.update(KProcID(4), 0.7);
    ASSERT_EQ(queue.getEventTime(KProcID(4)), 0.7);
    auto first2 = queue.getFirst();
    ASSERT_EQ(first2.first, 0.25);
    ASSERT_EQ(first2.second.data(), 5);

    queue.update(KProcID(1), 1.5);
    queue.update(KProcID(2), 1.75);
    queue.update(KProcID(3), 1.85);
    queue.update(KProcID(4), 1.2);
    queue.update(KProcID(5), 1.25);

    ASSERT_EQ(queue.getEventTime(KProcID(1)), 1.5);
    ASSERT_EQ(queue.getEventTime(KProcID(2)), 1.75);
    ASSERT_EQ(queue.getEventTime(KProcID(3)), 1.85);
    ASSERT_EQ(queue.getEventTime(KProcID(4)), 1.2);
    ASSERT_EQ(queue.getEventTime(KProcID(5)), 1.25);

    auto first3 = queue.getFirst();
    ASSERT_EQ(first3.first, std::numeric_limits<double>::infinity());
    ASSERT_EQ(first3.second.data(), 0);
}

TEST(KProcEventQueue, TimeCollision) {
    EventQueue queue;
    queue.updateMaxTime(1);
    queue.update(KProcID(2), 0.75);
    queue.update(KProcID(3), 0.75);
    queue.update(KProcID(4), 0.75);

    auto first1 = queue.getFirst();
    ASSERT_EQ(first1.first, 0.75);
    ASSERT_EQ(first1.second.data(), 2);

    queue.update(KProcID(2), 0.85);

    auto first2 = queue.getFirst();
    ASSERT_EQ(first2.first, 0.75);
    ASSERT_EQ(first2.second.data(), 3);

    queue.update(KProcID(3), 0.95);

    auto first3 = queue.getFirst();
    ASSERT_EQ(first3.first, 0.75);
    ASSERT_EQ(first3.second.data(), 4);
}

int main(int argc, char **argv) {
    int r = 0;
    ::testing::InitGoogleTest(&argc, argv);
    r = RUN_ALL_TESTS();
    return r;
}
