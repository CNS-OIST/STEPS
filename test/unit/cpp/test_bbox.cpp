#include <cmath>
#include <iterator>
#include <numeric>

#include "math/bbox.hpp"
#include "math/point.hpp"

#include <catch2/catch_test_macros.hpp>

using steps::math::bounding_box;
using steps::math::point3d;

TEST_CASE("BoundingBox_constructors") {
    bounding_box b_empty;
    REQUIRE(b_empty.empty());

    point3d pzero{0, 0, 0}, pones{1, 1, 1};
    REQUIRE_FALSE(b_empty.contains(pzero));

    bounding_box b0(pzero);
    REQUIRE_FALSE(b0.empty());
    REQUIRE(b0.min() == pzero);
    REQUIRE(b0.max() == pzero);

    bounding_box b00(pzero, pzero);
    REQUIRE_FALSE(b00.empty());
    REQUIRE(b00.min() == pzero);
    REQUIRE(b00.max() == pzero);

    bounding_box b01(pzero, pones);
    REQUIRE_FALSE(b0.empty());
    REQUIRE(b01.min() == pzero);
    REQUIRE(b01.max() == pones);

    bounding_box b10(pones, pzero);
    REQUIRE(b10.empty());

    point3d p[] = {{0, 0, 0}, {1, 0, 0.5}, {0, 0.5, 2}};
    bounding_box b3(std::begin(p), std::end(p));
    REQUIRE_FALSE(b3.empty());
    REQUIRE(b3.min() == pzero);
    REQUIRE(b3.max() == point3d(1, 0.5, 2));
}

TEST_CASE("BoundingBox_insert") {
    bounding_box b;

    point3d p[] = {{0, 0.2, 0}, {1, -0.5, 0.5}, {0, 0.5, 2}};

    for (const auto& x: p) {
        b.insert(x);
    }
    REQUIRE_FALSE(b.empty());

    REQUIRE(b.min() == point3d(0, -0.5, 0));
    REQUIRE(b.max() == point3d(1, 0.5, 2));
}

TEST_CASE("BoundingBox_contains") {
    constexpr size_t n = 3;
    point3d p[n] = {{0, 0.2, 0}, {1, -0.5, 0.5}, {0, 0.5, 2}};
    bounding_box b(p, p + n);

    REQUIRE_FALSE(b.empty());

    for (const auto& x: p) {
        REQUIRE(b.contains(x));
    }

    point3d pmean = std::accumulate(&p[0], &p[0] + n, point3d{0, 0, 0}) / static_cast<double>(n);

    REQUIRE(b.contains(pmean));

    REQUIRE(b.contains(b.max()));
    REQUIRE_FALSE(b.contains(b.max() + point3d(0.1, 0, 0)));
    REQUIRE_FALSE(b.contains(b.max() + point3d(0, 0.1, 0)));
    REQUIRE_FALSE(b.contains(b.max() + point3d(0, 0, 0.1)));

    REQUIRE(b.contains(b.min()));
    REQUIRE_FALSE(b.contains(b.min() - point3d(0.1, 0, 0)));
    REQUIRE_FALSE(b.contains(b.min() - point3d(0, 0.1, 0)));
    REQUIRE_FALSE(b.contains(b.min() - point3d(0, 0, 0.1)));
}

TEST_CASE("BoundingBox_clear") {
    point3d p[] = {{0, 0.2, 0}, {1, -0.5, 0.5}, {0, 0.5, 2}};
    bounding_box b(std::begin(p), std::end(p));

    b.clear();
    REQUIRE(b.empty());
    for (const auto& x: p) {
        REQUIRE_FALSE(b.contains(x));
    }

    b.insert(p[0]);
    REQUIRE(b.min() == p[0]);
    REQUIRE(b.max() == p[0]);
}
