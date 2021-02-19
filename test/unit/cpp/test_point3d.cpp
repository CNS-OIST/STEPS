#include <limits>
#include <cmath>

#include "steps/math/point.hpp"
#include "steps/math/tetrahedron.hpp"
#include "steps/math/triangle.hpp"

#include "gtest/gtest.h"

using namespace steps::math;

TEST(Point3d,constructors) {
    point3d p1{-1};

    ASSERT_EQ(p1[0],-1.0);

    point3d p2{-1,-2};

    ASSERT_EQ(p2[0],-1.0);
    ASSERT_EQ(p2[1],-2.0);

    point3d p3a{-1,-2,-3};
    point3d p3b={-1,-2,-3};

    ASSERT_EQ(p3a[0],-1.0);
    ASSERT_EQ(p3a[1],-2.0);
    ASSERT_EQ(p3a[2],-3.0);

    ASSERT_EQ(p3b[0],-1.0);
    ASSERT_EQ(p3b[1],-2.0);
    ASSERT_EQ(p3b[2],-3.0);
}

TEST(Point3d,indexed_access) {
    point3d p={2,4,8};

    std::swap(p[0],p[2]);
    p[1]*=3;

    ASSERT_EQ(p[0],8.0);
    ASSERT_EQ(p[1],12.0);
    ASSERT_EQ(p[2],2.0);
}

TEST(Point3d,arithmetic) {
    // expect exact arithmetic on 32-bit int values

    point3d p={2,4,8};

    point3d q={1,3,5};
    point3d q0=q;
    q+=p;
    ASSERT_EQ(q,point3d(3,7,13));

    q-=p;
    ASSERT_EQ(q,q0);

    ASSERT_EQ(point3d(1,2,3)+point3d(4,5,6),point3d(5,7,9));
    ASSERT_EQ(point3d(1,2,3)-point3d(6,5,4),point3d(-5,-3,-1));
}

TEST(Point3d,iterators) {
    point3d p={7,6,5};

    std::vector<double> v(p.begin(), p.end());
    ASSERT_EQ(p[0],v[0]);
    ASSERT_EQ(p[1],v[1]);
    ASSERT_EQ(p[2],v[2]);
}

TEST(Point3d,vector_ops) {
    point3d p={7,6,5};
    point3d q={0.5,0.25,0.125};

    double d=p[0]*q[0]+p[1]*q[1]+p[2]*q[2];

    ASSERT_EQ(d,p.dot(q));

    point3d x={p[1]*q[2]-p[2]*q[1],
               p[2]*q[0]-p[0]*q[2],
               p[0]*q[1]-p[1]*q[0]};

    ASSERT_EQ(x,p.cross(q));
}

TEST(Point3d,normL2) {
    point3d p={1,-2,3};

    auto expected = sqrt(p[0]*p[0]+p[1]*p[1]+p[2]*p[2]);
    ASSERT_DOUBLE_EQ(norm(p), expected);
}

TEST(Point3d,normL1) {
    point3d p={1,2,-3};

    auto expected = std::abs(p[0])+std::abs(p[1])+std::abs(p[2]);
    ASSERT_DOUBLE_EQ(normL1(p), expected);

    ASSERT_GE(normL1(p), norm(p));
}

TEST(Tetrahedron,tet_inside) {

    point3d t0{0, 0, 0},
            t1{2.2, 0, 0},
            t2{0, 3.3, 0},
            t3{0, 0, 4.4};

    // clearly inside
    {
        point3d p0{0.1, 0.1, 0.1};
        ASSERT_TRUE(tet_inside(t0, t1, t2, t3, p0));
    }
    // clearly outside
    {
        point3d p0{3.1, 0.1, 0.1};
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
    {
        point3d p0{0.1, 10.1, 0.1};
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
    {
        point3d p0{0.1, 0.1, 5.1};
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }

}

TEST(Tetrahedron,tet_inside_face_sampling) {

    // this fails reliably without tolerance in sameside

    point3d t0{0.0, 0, 0},
            t1{2.2, 0, 0},
            t2{0, 3.3, 0},
            t3{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100,200,300};

    t0 = t0*scale + move;
    t1 = t1*scale + move;
    t2 = t2*scale + move;
    t3 = t3*scale + move;

    {
        point3d e1 = t1-t3;
        point3d e2 = t2-t3;

        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            double alpha = static_cast<double>(n)/(nsample-1);
            for (size_t m = 0; m < nsample; ++m) {
                double beta = static_cast<double>(m)/(nsample-1);
                if (alpha+beta <= 1.0){
                    point3d p0 = e1*alpha + e2*beta + t3;
                    ASSERT_TRUE(tet_inside(t0, t1, t2, t3, p0));
                }
            }
        }
        // clearly outside
        point3d p0 = t1*0.3 + t2*0.5 + t3*2.0;
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }
}

TEST(Tetrahedron,tet_inside_edge_sampling) {

    // this fails reliably without tolerance in sameside

    point3d t0{1, 0, 0},
            t1{2.2, 0, 0},
            t2{0, 3.3, 0},
            t3{0, 0, 4.4};

    std::array<point3d,3> vertices{t1, t2, t3};

    for (const auto& v : vertices)
    {
        point3d diff = v-t0;
        size_t nsample = 100;
        for (size_t n = nsample/2; n < nsample; ++n) {
            point3d p0 = t0 + diff * n/(nsample-1);
            ASSERT_TRUE(tet_inside(t0, t1, t2, t3, p0));
        }
        // clearly outside
        point3d p0 = t0 + 1.3*diff;
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p0));
    }

}

TEST(Tetrahedron,tet_inside_vertices) {

    point3d t0{-1, 0, 0},
            t1{2.2, 0, 0},
            t2{0, 3.3, 0},
            t3{0, 0, 4.4};

    std::array<point3d,4> vertices{t0, t1, t2, t3};

    for (const auto& v : vertices)
    {
        ASSERT_TRUE(tet_inside(t0, t1, t2, t3, v));
    }

}

TEST(Triangle,tri_normal) {

    point3d t0{0, 0, 0},
            t1{1, 0, 0},
            t2{0, 1, 0};

    point3d expected{0, 0, 1};
    point3d normal = tri_normal(t0, t1, t2);

    ASSERT_DOUBLE_EQ(normal[0], expected[0]);
    ASSERT_DOUBLE_EQ(normal[1], expected[1]);
    ASSERT_DOUBLE_EQ(normal[2], expected[2]);
}

TEST(Triangle,parallel2face) {

    point3d t0{2.2, 0, 0},
            t1{0, 3.3, 0},
            t2{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100,200,300};

    t0 = t0*scale + move;
    t1 = t1*scale + move;
    t2 = t2*scale + move;

    point3d lc(tri_barycenter(t0,t1,t2));

    std::array<point3d,2> vertices{t1, t2};

    point3d intersect;

    for (const auto& v : vertices)
    {
        // point between t0 and v
        point3d pm = t0 + 0.3*(v-t0);
        // edge between barycenter and midpoint
        point3d diff = pm-lc;
        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            // this point is outside
            point3d p1 = lc + diff * (1.5+static_cast<double>(n)/(nsample-1));
            // p1-lc is in face plane but not parallel to any edge
            ASSERT_FALSE(tri_intersect_line(t0, t1, t2, lc, p1, intersect));
        }
    }
}

TEST(Triangle,intersection_on_edge) {

    point3d t0{2.2, 0, 0},
            t1{0, 3.3, 0},
            t2{0, 0, 4.4};

    // scale
    const double scale = 50;
    // translate
    const point3d move{100,200,300};

    t0 = t0*scale + move;
    t1 = t1*scale + move;
    t2 = t2*scale + move;

    std::array<point3d,2> vertices{t1, t2};

    point3d p0{0,0,0};
    point3d intersect;

    for (const auto& v : vertices)
    {
        point3d diff = v-t0;
        size_t nsample = 100;
        for (size_t n = 0; n < nsample; ++n) {
            point3d p1 = t0 + diff * n/(nsample-1);
            ASSERT_TRUE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
        }
        // endpoint is outside the edge
        point3d p1 = t0 + 1.3*diff;
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
    }
}

TEST(Triangle,intersection) {

    point3d intersect;
    {
        point3d t0{1, 0, 0},
                t1{0, 1, 0},
                t2{0, 0, 1};

        point3d lc(tri_barycenter(t0,t1,t2)),
                l1{100.0, 100.0, 100.0},
                l0(lc - 0.2*(l1-lc));

        // Common case
        ASSERT_TRUE(tri_intersect_line(t0, t1, t2, l0, l1, intersect));
        ASSERT_TRUE(tri_intersect_line(t0, t1, t2, l0, l1, intersect, false));
        ASSERT_TRUE(lc.almostEqual(intersect));

        // outside yields false
        point3d l2 = lc + 1.5*(l1-lc);  // point after l1
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, l1, l2, intersect, true));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, l1, l2, intersect, false));

        //But reversal ok in ray mode
        ASSERT_TRUE(tri_intersect_line(t0, t1, t2, l2, l1, intersect, false));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, l2, l1, intersect, true));

        {
            point3d p0{0.0, 2.0, 1.0}, p1{2.0, 2.0, 1.0};
            ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
            ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
        }
        {
            // 100% Independent
            point3d p0{0.0, 0.2, 0.2},
                    p1{0.2, 0.2, 0.2};
            ASSERT_TRUE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
            ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
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
        ASSERT_TRUE(tet_inside(t0, t1, t2, t3, p0));
        // Check that p1 is outside the tet (t0, t1, t2, t3)
        ASSERT_FALSE(tet_inside(t0, t1, t2, t3, p1));
        // So [p0, p1] should intersect the triangle (t1, t2, t3)
        ASSERT_TRUE(tri_intersect_line(t1, t2, t3, p0, p1, intersect));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t3, p0, p1, intersect));
        ASSERT_FALSE(tri_intersect_line(t0, t2, t3, p0, p1, intersect));
    }
}
