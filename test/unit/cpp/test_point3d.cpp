#include <limits>
#include <cmath>

#include "steps/math/point.hpp"
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

TEST(Point3d,intersection) {
    point3d t0{1, 1, 1},
            t1{1, 1, 2},
            t2{1, 2, 2};

    // result holder
    point3d intersect;

    {
        point3d lc(tri_barycenter(t0,t1,t2)),
                l1{10.0, 10.0, 10.0},
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

        point3d p0{0.0, 2.0, 1.0}, p1{2.0, 2.0, 1.0};
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
    }

    // 100% Independent
    {
        point3d p0{0.0, 1.2, 1.8},
                p1{0.5, 1.2, 1.8};
        ASSERT_TRUE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, false));
        ASSERT_FALSE(tri_intersect_line(t0, t1, t2, p0, p1, intersect, true));
    }

}
