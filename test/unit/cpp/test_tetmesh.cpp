#include <iostream>
#include <memory>
#include <limits>
#include <cmath>

#include "steps/geom/tetmesh.hpp"

#include "gtest/gtest.h"

#define COORDS v_coords
#define TETINDICES t_indices
#include "./sample_meshdata.h"
#undef COORDS
#undef TETINDICES

using steps::index_t;
using steps::tetmesh::Tetmesh;
using steps::math::point3d;

struct TetmeshTest: public ::testing::Test {
    std::unique_ptr<Tetmesh> mesh;
    size_t vsN;
    size_t tN;

    virtual void SetUp() {
        const double *vs=&v_coords[0][0];
        vsN=sizeof(v_coords)/sizeof(*vs);

        const auto *ts=&t_indices[0][0];
        tN=sizeof(t_indices)/sizeof(*ts);

        mesh.reset(new Tetmesh(std::vector<double>(vs, vs + vsN),
                               std::vector<steps::vertex_id_t::value_type>(ts, ts + tN)));
    }
};

TEST_F(TetmeshTest,user_constructor) {
    ASSERT_EQ(mesh->countVertices(),vsN/3);
    ASSERT_EQ(mesh->countTets(),tN/4);
}

TEST_F(TetmeshTest,tri_norm) {
    size_t triN = mesh->countTris();

    for (index_t i = 0u; i < triN; ++i) {
        const auto *tri = mesh->_getTri(i);

        point3d v[]={mesh->_getVertex(tri[0]),
                     mesh->_getVertex(tri[1]),
                     mesh->_getVertex(tri[2])};

        point3d n=(v[1]-v[0]).cross(v[2]-v[0]);
        n /= std::sqrt(n.dot(n));

        auto m = mesh->getTriNorm(i);
        ASSERT_DOUBLE_EQ(n[0],m[0]);
        ASSERT_DOUBLE_EQ(n[1],m[1]);
        ASSERT_DOUBLE_EQ(n[2],m[2]);
    }
}



TEST_F(TetmeshTest,intersect) {
    double v[] = {0.5, 0.5, 0., 0.5, 0.5, 1.5};
    point3d p0(v[0], v[1], v[2]),
            p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    ASSERT_FALSE(res.empty());
    // Intersects tet 0,1,2
    ASSERT_EQ(res[0].first, index_t(0));
    ASSERT_GT(res[0].second, 0.04);

    // Batch mode
    auto res2 = mesh->intersect(v, 2);
    for (const auto& t: res2)
        for(const auto& r: t)
            std::cout << r.first << ": " << r.second << std::endl;
    ASSERT_EQ(res2.size(), 1ul);
    ASSERT_EQ(res2[0].size(), 3ul);
    ASSERT_EQ(res2[0][0].first, 0ul);
    ASSERT_GT(res2[0][0].second, 0.04);

}


TEST_F(TetmeshTest,intersect2) {
    Tetmesh tm({0.,0.,0.,5.,5.,0.,5.,0.,0.,0.,0.,5.},
               {0,1,2,3});

    auto res = tm.intersectDeterministic(point3d(.0, .0, .0), point3d(0.01, 0.01, 4.));
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Intersects tet 0 for almost 100%
    ASSERT_EQ(res.size(), 1ul);
    ASSERT_EQ(res[0].first, index_t(0));
    ASSERT_GT(res[0].second, 0.90);

}
