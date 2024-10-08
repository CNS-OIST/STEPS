#include "util/pqueue.hpp"

#include <catch2/catch_test_macros.hpp>


TEST_CASE("pqueue") {
    steps::util::PQueue<double, int> pq;
    // add in non-sorted order
    pq.push(0, 0);
    pq.push(3.1, 3);
    pq.push(2, 1);
    pq.push(2.5, 2);
    // add twice
    pq.push(4, 1);
    // override processed of element in the queue
    pq.mark_as_processed(2);
    // override processed of element not in the queue
    pq.mark_as_processed(5);
    pq.push(5, 5);

    std::vector<int> res = {0, 1, 3};

    size_t ii = 0;
    while (!pq.empty()) {
        REQUIRE(pq.next().second == res[ii]);
        // try to add at the end of the queue.
        pq.push(5, 1);
        ++ii;
    }
    REQUIRE(ii == res.size());
    REQUIRE(pq.get_min_unqueued_value() == 4);
    pq.push(4, 4);
    REQUIRE(pq.get_min_unqueued_value() == 6);
    REQUIRE(pq.size_processed() == 5);
}
