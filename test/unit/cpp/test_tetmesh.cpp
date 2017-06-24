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

using steps::tetmesh::Tetmesh;
using steps::math::point3d;

struct TetmeshTest: public ::testing::Test {
    std::unique_ptr<Tetmesh> mesh;
    size_t vsN;
    size_t tN;

    virtual void SetUp() {
        const double *vs=&v_coords[0][0];
        vsN=sizeof(v_coords)/sizeof(*vs);

        const unsigned int *ts=&t_indices[0][0];
        tN=sizeof(t_indices)/sizeof(*ts);

        mesh.reset(new Tetmesh(std::vector<double>(vs,vs+vsN),std::vector<unsigned int>(ts,ts+tN)));
    }
};

TEST_F(TetmeshTest,user_constructor) {
    ASSERT_EQ(mesh->countVertices(),vsN/3);
    ASSERT_EQ(mesh->countTets(),tN/4);
}

TEST_F(TetmeshTest,tri_norm) {
    size_t triN = mesh->countTris();

    for (size_t i=0; i<triN; ++i) {
        const uint *tri = mesh->_getTri(i);

        point3d v[]={mesh->_getVertex(tri[0]),
                     mesh->_getVertex(tri[1]),
                     mesh->_getVertex(tri[2])};

        point3d n=cross(v[1]-v[0],v[2]-v[0]);
        n /= std::sqrt(dot(n,n));

        auto m = mesh->getTriNorm(i);
        ASSERT_DOUBLE_EQ(n[0],m[0]);
        ASSERT_DOUBLE_EQ(n[1],m[1]);
        ASSERT_DOUBLE_EQ(n[2],m[2]);
    }
}
