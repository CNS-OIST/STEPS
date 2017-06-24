#include <iterator>
#include <numeric>
#include <cmath>

#include "steps/math/point.hpp"
#include "steps/math/bbox.hpp"

#include "gtest/gtest.h"

using namespace steps::math;

TEST(BoundingBox,constructors) {
    bounding_box b_empty;
    ASSERT_TRUE(b_empty.empty());

    point3d pzero{0,0,0}, pones{1,1,1};
    ASSERT_FALSE(b_empty.contains(pzero));

    bounding_box b0(pzero);
    ASSERT_FALSE(b0.empty());
    ASSERT_EQ(b0.min(),pzero);
    ASSERT_EQ(b0.max(),pzero);

    bounding_box b00(pzero,pzero);
    ASSERT_FALSE(b00.empty());
    ASSERT_EQ(b00.min(),pzero);
    ASSERT_EQ(b00.max(),pzero);

    bounding_box b01(pzero,pones);
    ASSERT_FALSE(b0.empty());
    ASSERT_EQ(b01.min(),pzero);
    ASSERT_EQ(b01.max(),pones);

    bounding_box b10(pones,pzero);
    ASSERT_TRUE(b10.empty());

    point3d p[]={{0,0,0},{1,0,0.5},{0,0.5,2}};
    bounding_box b3(std::begin(p),std::end(p));
    ASSERT_FALSE(b3.empty());
    ASSERT_EQ(b3.min(),pzero);
    ASSERT_EQ(b3.max(),point3d(1,0.5,2));
}

TEST(BoundingBox,insert) {
    bounding_box b;

    point3d p[]={{0,0.2,0},{1,-0.5,0.5},{0,0.5,2}};

    for (const auto &x: p) b.insert(x);
    ASSERT_FALSE(b.empty());

    ASSERT_EQ(b.min(),point3d(0,-0.5,0));
    ASSERT_EQ(b.max(),point3d(1,0.5,2));
}

TEST(BoundingBox,contains) {
    constexpr size_t n=3;
    point3d p[n]={{0,0.2,0},{1,-0.5,0.5},{0,0.5,2}};
    bounding_box b(p,p+n);

    ASSERT_FALSE(b.empty());

    for (const auto &x: p) {
        ASSERT_TRUE(b.contains(x));
    }

    point3d pmean=std::accumulate(&p[0],&p[0]+n,point3d{0,0,0})/(double)n;

    ASSERT_TRUE(b.contains(pmean));

    ASSERT_TRUE(b.contains(b.max()));
    ASSERT_FALSE(b.contains(b.max()+point3d(0.1,0,0)));
    ASSERT_FALSE(b.contains(b.max()+point3d(0,0.1,0)));
    ASSERT_FALSE(b.contains(b.max()+point3d(0,0,0.1)));

    ASSERT_TRUE(b.contains(b.min()));
    ASSERT_FALSE(b.contains(b.min()-point3d(0.1,0,0)));
    ASSERT_FALSE(b.contains(b.min()-point3d(0,0.1,0)));
    ASSERT_FALSE(b.contains(b.min()-point3d(0,0,0.1)));
}

TEST(BoundingBox,clear) {
    point3d p[]={{0,0.2,0},{1,-0.5,0.5},{0,0.5,2}};
    bounding_box b(std::begin(p),std::end(p));

    b.clear();
    ASSERT_TRUE(b.empty());
    for (const auto &x: p) {
        ASSERT_FALSE(b.contains(x));
    }

    b.insert(p[0]);
    ASSERT_EQ(b.min(),p[0]);
    ASSERT_EQ(b.max(),p[0]);
}
