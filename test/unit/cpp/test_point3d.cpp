#include <cmath>
#include <limits>

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "math/point.hpp"
#include "math/segment.hpp"
#include "math/tetrahedron.hpp"
#include "math/triangle.hpp"

using namespace steps::math;
// these segments are in the order of 1 so this is also the linear tolerance
constexpr double adim_tol = 4 * std::numeric_limits<point3d::value_type>::epsilon();
constexpr double eps = adim_tol / 2.0;

TEST_CASE("Point3d_constructors") {
    point3d p1;

    REQUIRE(p1[0] == 0.0);
    REQUIRE(p1[1] == 0.0);
    REQUIRE(p1[2] == 0.0);

    point3d p3a{-1, -2, -3};
    point3d p3b = {-1, -2, -3};

    REQUIRE(p3a[0] == -1.0);
    REQUIRE(p3a[1] == -2.0);
    REQUIRE(p3a[2] == -3.0);

    REQUIRE(p3b[0] == -1.0);
    REQUIRE(p3b[1] == -2.0);
    REQUIRE(p3b[2] == -3.0);
}

TEST_CASE("Point3d_indexed_access") {
    point3d p = {2, 4, 8};

    std::swap(p[0], p[2]);
    p[1] *= 3;

    REQUIRE(p[0] == 8.0);
    REQUIRE(p[1] == 12.0);
    REQUIRE(p[2] == 2.0);
}

TEST_CASE("Point3d_arithmetic") {
    // expect exact arithmetic on 32-bit int values

    point3d p = {2, 4, 8};

    point3d q = {1, 3, 5};
    point3d q0 = q;
    q += p;
    REQUIRE(q == point3d(3, 7, 13));

    q -= p;
    REQUIRE(q == q0);

    REQUIRE(point3d(1, 2, 3) + point3d(4, 5, 6) == point3d(5, 7, 9));
    REQUIRE(point3d(1, 2, 3) - point3d(6, 5, 4) == point3d(-5, -3, -1));
}

TEST_CASE("Point3d_iterators") {
    point3d p = {7, 6, 5};

    std::vector<double> v(p.begin(), p.end());
    REQUIRE(p[0] == v[0]);
    REQUIRE(p[1] == v[1]);
    REQUIRE(p[2] == v[2]);
}

TEST_CASE("Point3d_vector_ops") {
    point3d p = {7, 6, 5};
    point3d q = {0.5, 0.25, 0.125};

    double d = p[0] * q[0] + p[1] * q[1] + p[2] * q[2];

    REQUIRE(d == p.dot(q));

    point3d x = {p[1] * q[2] - p[2] * q[1], p[2] * q[0] - p[0] * q[2], p[0] * q[1] - p[1] * q[0]};

    REQUIRE(x == p.cross(q));
}

TEST_CASE("Point3d_normL2") {
    point3d p = {1, -2, 3};

    auto expected = sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
    REQUIRE_THAT(norm(p), Catch::Matchers::WithinULP(expected, 4));
}

TEST_CASE("Point3d_normL1") {
    point3d p = {1, 2, -3};

    auto expected = std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2]);
    REQUIRE_THAT(normL1(p), Catch::Matchers::WithinULP(expected, 4));

    REQUIRE(normL1(p) >= norm(p));
}

TEST_CASE("Tetrahedron_tet_inside") {
    point3d t0{0, 0, 0}, t1{2.2, 0, 0}, t2{0, 3.3, 0}, t3{0, 0, 4.4};

    // clearly inside
    {
        point3d p0{0.1, 0.1, 0.1};
        REQUIRE(tet_inside(t0, t1, t2, t3, p0));
    }
    // clearly outside
    {
        point3d p0{3.1, 0.1, 0.1};
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
    {
        point3d p0{0.1, 10.1, 0.1};
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
    {
        point3d p0{0.1, 0.1, 5.1};
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
}

TEST_CASE("Tetrahedron_tet_inside_face_sampling") {
    // this fails reliably without tolerance in sameside

    point3d t0{0.0, 0, 0}, t1{2.2, 0, 0}, t2{0, 3.3, 0}, t3{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100, 200, 300};

    t0 = t0 * scale + move;
    t1 = t1 * scale + move;
    t2 = t2 * scale + move;
    t3 = t3 * scale + move;

    {
        point3d e1 = t1 - t3;
        point3d e2 = t2 - t3;

        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            double alpha = static_cast<double>(n) / (nsample - 1);
            for (size_t m = 0; m < nsample; ++m) {
                double beta = static_cast<double>(m) / (nsample - 1);
                if (alpha + beta <= 1.0) {
                    point3d p0 = e1 * alpha + e2 * beta + t3;
                    REQUIRE(tet_inside(t0, t1, t2, t3, p0));
                }
            }
        }
        // clearly outside
        point3d p0 = t1 * 0.3 + t2 * 0.5 + t3 * 2.0;
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
}

TEST_CASE("Tetrahedron_tet_inside_edge_sampling") {
    // this fails reliably without tolerance in sameside

    point3d t0{1, 0, 0}, t1{2.2, 0, 0}, t2{0, 3.3, 0}, t3{0, 0, 4.4};

    std::array<point3d, 3> vertices{t1, t2, t3};

    for (const auto& v: vertices) {
        point3d diff = v - t0;
        size_t nsample = 100;
        for (size_t n = nsample / 2; n < nsample; ++n) {
            point3d p0 = t0 + diff * n / (nsample - 1);
            REQUIRE(tet_inside(t0, t1, t2, t3, p0));
        }
        // clearly outside
        point3d p0 = t0 + 1.3 * diff;
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
}

TEST_CASE("Tetrahedron_tet_inside_vertices") {
    point3d t0{-1, 0, 0}, t1{2.2, 0, 0}, t2{0, 3.3, 0}, t3{0, 0, 4.4};

    std::array<point3d, 4> vertices{t0, t1, t2, t3};

    for (const auto& v: vertices) {
        REQUIRE(tet_inside(t0, t1, t2, t3, v));
    }
}

TEST_CASE("Triangle_tri_normal") {
    point3d t0{0, 0, 0}, t1{1, 0, 0}, t2{0, 1, 0};

    point3d expected{0, 0, 1};
    point3d normal = tri_normal(t0, t1, t2);

    REQUIRE_THAT(normal[0], Catch::Matchers::WithinULP(expected[0], 4));
    REQUIRE_THAT(normal[1], Catch::Matchers::WithinULP(expected[1], 4));
    REQUIRE_THAT(normal[2], Catch::Matchers::WithinULP(expected[2], 4));
}

TEST_CASE("Triangle_parallel2face") {
    point3d t0{2.2, 0, 0}, t1{0, 3.3, 0}, t2{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100, 200, 300};

    t0 = t0 * scale + move;
    t1 = t1 * scale + move;
    t2 = t2 * scale + move;

    point3d lc(tri_barycenter(t0, t1, t2));

    std::array<point3d, 2> vertices{t1, t2};

    point3d intersect;

    for (const auto& v: vertices) {
        // point between t0 and v
        point3d pm = t0 + 0.3 * (v - t0);
        // edge between barycenter and midpoint
        point3d diff = pm - lc;
        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            // this point is outside
            point3d p1 = lc + diff * (1.5 + static_cast<double>(n) / (nsample - 1));
            // p1-lc is in face plane but not parallel to any edge
            REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, lc, p1, intersect));
        }
    }
}

TEST_CASE("Triangle_intersection_on_edge") {
    point3d t0{2.2, 0, 0}, t1{0, 3.3, 0}, t2{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100, 200, 300};

    t0 = t0 * scale + move;
    t1 = t1 * scale + move;
    t2 = t2 * scale + move;

    std::array<point3d, 2> vertices{t1, t2};

    point3d p0{0, 0, 0};
    point3d intersect;

    for (const auto& v: vertices) {
        point3d diff = v - t0;
        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            point3d p1 = t0 + diff * n / (nsample - 1);
            REQUIRE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
        }
        // endpoint is outside the edge
        point3d p1 = t0 + 1.3 * diff;
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
    }
}

TEST_CASE("Triangle_intersection_line") {
    point3d intersect;
    {
        point3d t0{1, 0, 0}, t1{0, 1, 0}, t2{0, 0, 1};
        point3d lc(tri_barycenter(t0, t1, t2)), l1{100.0, 100.0, 100.0}, l0(lc - 0.2 * (l1 - lc));

        // Common case
        REQUIRE(tri_intersect_line(t0, t1, t2, l0, l1, intersect));
        REQUIRE(tri_intersect_line(t0, t1, t2, l0, l1, intersect, false));
        REQUIRE(lc.almostEqual(intersect, adim_tol * l0.distance(l1)));

        // outside yields false
        point3d l2 = lc + 1.5 * (l1 - lc);  // point after l1
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, l1, l2, intersect, true));
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, l1, l2, intersect, false));

        // But reversal ok in ray mode
        REQUIRE(tri_intersect_line(t0, t1, t2, l2, l1, intersect, false));
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, l2, l1, intersect, true));

        {
            point3d p0{0.0, 2.0, 1.0}, p1{2.0, 2.0, 1.0};
            REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
            REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
        }
        {
            // 100% Independent
            point3d p0{0.0, 0.2, 0.2}, p1{0.2, 0.2, 0.2};
            REQUIRE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
            REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
        }
    }

    // p0 is on the triangle (t1, t2, t3), p1 outside
    {
        point3d p0{-41.965698240000002, 155.23899840999999, 23.905799869999999};
        point3d p1{-46.6753006, 165.61300659, 28.899000170000001};

        point3d t0{52.625232766555143, -181.95954873227998, -2.5301888029974244};
        point3d t1{-28.186250442902992, 172.50887046333756, 37.827370072762626};
        point3d t2{-92.357735303495843, 89.356897035748702, -308.63369376714121};
        point3d t3{-202.00328266203834, -43.904510678466544, -67.190066180242553};

        // Check that p0 is inside the tet (t0, t1, t2, t3)
        REQUIRE(tet_inside(t0, t1, t2, t3, p0));
        // Check that p1 is outside the tet (t0, t1, t2, t3)
        REQUIRE_FALSE(tet_inside(t0, t1, t2, t3, p1));
        // So [p0, p1] should intersect the triangle (t1, t2, t3)
        REQUIRE(tri_intersect_line(t1, t2, t3, p0, p1, intersect));
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
        REQUIRE_FALSE(tri_intersect_line(t0, t1, t3, p0, p1, intersect));
        REQUIRE_FALSE(tri_intersect_line(t0, t2, t3, p0, p1, intersect));
    }
}

TEST_CASE("Triangle_intersection_point") {
    point3d t0{0, 0, 0}, t1{0, 1, 0}, t2{1, 0, 0};
    const auto eps2 = std::numeric_limits<point3d::value_type>::epsilon() / 2.0;
    // inside
    REQUIRE(tri_intersect_point(t0, t1, t2, {0.1, 0.1, 0}));
    // inversion
    REQUIRE(tri_intersect_point(t1, t0, t2, {0.1, 0.1, 0}));
    // edge
    REQUIRE(tri_intersect_point(t0, t1, t2, {0.1, 0, 0}));
    // vertex
    REQUIRE(tri_intersect_point(t0, t1, t2, {0, 0, 0}));
    // vertex+eps (inside)
    REQUIRE(tri_intersect_point(t0, t1, t2, {eps2, eps2, eps2}));
    // vertex-eps (outside but it should be inside)
    REQUIRE(tri_intersect_point(t0, t1, t2, {-eps2, -eps2, -eps2}));
    // edge+eps (inside)
    REQUIRE(tri_intersect_point(t0, t1, t2, {0.1, eps2, 0}));
    // edge-eps (outside but it should be inside)
    REQUIRE(tri_intersect_point(t0, t1, t2, {0.1, -eps2, 0}));
    // outside 2D
    REQUIRE(!tri_intersect_point(t0, t1, t2, {1, 1, 0}));
    // outside 3D
    REQUIRE(!tri_intersect_point(t0, t1, t2, {0, 0, 1}));
    // inside 3D+eps (outside but it should be inside)
    REQUIRE(tri_intersect_point(t0, t1, t2, {0.1, 0.1, eps2}));
}

TEST_CASE("Segment_point_intersection") {
    point3d a{0, 0, 0}, b{1, 0, 0};
    REQUIRE(segment_intersect_point(a, b, {0.5, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {0.5 + eps, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {0.5 + eps, eps, eps}));
    REQUIRE(segment_intersect_point(a, b, {0.7, 0, eps}));
    REQUIRE(segment_intersect_point(a, b, {1.0 - eps, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {1.0, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {0, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {-eps, 0, 0}));
    REQUIRE(segment_intersect_point(a, b, {1.0 + eps, 0, 0}));
    REQUIRE(!segment_intersect_point(a, b, {-3.0 * eps, 0, 0}));
    REQUIRE(!segment_intersect_point(a, b, {0.7, 0, 3.0 * eps}));
    REQUIRE(segment_intersect_point(a, b, {-3.0 * eps, 0, 0}, false));
    REQUIRE(!segment_intersect_point(a, b, {0.7, 0, 3.0 * eps}, false));
}

TEST_CASE("Segment_segment_intersection") {
    point3d p0{0, 0, 0}, p1{1, 0, 0}, q0{0.25, 0.25, 0};
    const point3d reset_point{-1, -1, -1};

    // perfect intersection
    point3d intersection0 = reset_point;
    point3d intersection1 = reset_point;
    auto ans =
        segment_intersect_segment(p0, p1, q0, {0.25, -0.25, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({0.25, 0, 0}, adim_tol)));
    // skewed
    intersection0 = reset_point;
    ans = segment_intersect_segment(p0, p1, q0, {0.75, -0.25, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({0.5, 0, 0}, adim_tol)));
    // skewed, not coplanar but in tolerance
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {0.25, 0.25, eps}, {0.75, -0.25, eps}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({0.5, 0, 0}, adim_tol)));
    // not coplanar
    intersection0 = reset_point;
    ans = segment_intersect_segment(p0, p1, q0, {0.75, -0.25, 1}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // not coplanar 2
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {0.75, 0.25, 1}, {0.75, -0.25, 1}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // too short q
    intersection0 = reset_point;
    ans = segment_intersect_segment(p0, p1, q0, {0.25, 0.1, 0}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // too short q 2
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {0.25, -0.1, 0}, {0.25, -0.5, 0}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // too short p
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {-0.1, 0.25, 0}, {-0.1, -0.25, 0}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // too short p 2
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {1, 0.25, 0}, {1.2, -0.25, 0}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // extreme connection
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {0.9, 0.25, 0}, {1.1, -0.25, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({1, 0, 0}, adim_tol)));
    // extreme connection 2
    intersection0 = reset_point;
    ans =
        segment_intersect_segment(p0, p1, {0.9, 0.25, 0}, {1, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({1, 0, 0}, adim_tol)));
    // parallel
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {0, 0, 4 * eps}, {1, 0, 4 * eps}, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // colinear 2
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0, 0, 0}, {0.5, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.5, 0, 0}, adim_tol)));
    // colinear 3
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0.1, 0, 0}, {0.5, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0.1, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.5, 0, 0}, adim_tol)));
    // colinear 4
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {-2, 0, 0}, {-1, 0, 0}, intersection0, intersection1);
    REQUIRE(!ans);
    // colinear 5
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {3, 0, 0}, {2, 0, 0}, intersection0, intersection1);
    REQUIRE(!ans);
    // colinear 6
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {-1, 0, 0}, {0.2, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.2, 0, 0}, adim_tol)));
    // colinear 7
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0.2, 0, 0}, {-1, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.2, 0, 0}, adim_tol)));
    // colinear 8
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {-1, 0, 0}, {0, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0, 0, 0}, adim_tol)));
    // colinear 9
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {2, 0, 0}, {1, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{1, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{1, 0, 0}, adim_tol)));
    // colinear 10
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0.5, 0, 0}, {0.2, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0.2, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.5, 0, 0}, adim_tol)));
    // colinear 11
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {-1, 0, 0}, {2, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{1, 0, 0}, adim_tol)));
    // colinear 12
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {2, 0, 0}, {-1, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{1, 0, 0}, adim_tol)));
    // colinear 13
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {1, 0, 0}, {0, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{1, 0, 0}, adim_tol)));
    // colinear 14
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0.8, 0, 0}, {0, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual(point3d{0, 0, 0}, adim_tol)));
    REQUIRE((ans && intersection1.almostEqual(point3d{0.8, 0, 0}, adim_tol)));
    // not too short q if it is a line
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, q0, {0.25, 0.1, 0}, intersection0, intersection1, true, false);
    REQUIRE((ans == 1 && intersection0.almostEqual({0.25, 0, 0}, adim_tol)));
    // not too short p if it is a line
    intersection0 = reset_point;
    ans = segment_intersect_segment(
        p0, p1, {-0.1, 0.25, 0}, {-0.1, -0.25, 0}, intersection0, intersection1, false, true);
    REQUIRE((ans == 1 && intersection0.almostEqual({-0.1, 0, 0}, adim_tol)));
    // q has no length and does not lie on p
    intersection0 = reset_point;
    ans = segment_intersect_segment(p0, p1, q0, q0, intersection0, intersection1);
    REQUIRE((!ans && intersection0.almostEqual(reset_point, adim_tol)));
    // q has no length and lies on p
    intersection0 = reset_point;
    ans = segment_intersect_segment(p0, p1, {0.5, 0, 0}, {0.5, 0, 0}, intersection0, intersection1);
    REQUIRE((ans && intersection0.almostEqual({0.5, 0, 0}, adim_tol)));

    // skewed. Real life error
    intersection0 = reset_point;
    intersection1 = reset_point;
    ans = segment_intersect_segment({0.00012993, 0.00015986, 0.00022993},
                                    {0.0002, 0.0003, 0.0003},
                                    {9.99e-05, 0.0003, 0.0001999},
                                    {9.99e-05, 0.0001999, 0.0003},
                                    intersection0,
                                    intersection1);
    REQUIRE(!ans);
}
