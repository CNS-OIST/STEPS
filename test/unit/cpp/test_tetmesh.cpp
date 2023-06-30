#include "geom/tetmesh.hpp"

#include <cmath>
#include <iostream>
#include <limits>
#include <memory>

#include "gtest/gtest.h"

#define COORDS     v_coords
#define TETINDICES t_indices
#include "sample_meshdata.h"
#undef COORDS
#undef TETINDICES

using steps::index_t;
using steps::tetrahedron_global_id;
using steps::triangle_global_id;
using steps::math::point3d;
using steps::tetmesh::Tetmesh;

struct TetmeshTest: public ::testing::Test {
    std::unique_ptr<Tetmesh> mesh;
    size_t vsN;
    size_t tN;

    virtual void SetUp() {
        const double* vs = &v_coords[0][0];
        vsN = sizeof(v_coords) / sizeof(*vs);

        const auto* ts = &t_indices[0][0];
        tN = sizeof(t_indices) / sizeof(*ts);

        mesh.reset(new Tetmesh(std::vector<double>(vs, vs + vsN),
                               std::vector<steps::vertex_id_t::value_type>(ts, ts + tN)));
    }
};

TEST_F(TetmeshTest, user_constructor) {
    ASSERT_EQ(mesh->countVertices(), vsN / 3);
    ASSERT_EQ(mesh->countTets(), tN / 4);
}

TEST_F(TetmeshTest, tri_norm) {
    size_t triN = mesh->countTris();

    for (auto i: triangle_global_id::range(triN)) {
        const auto& tri = mesh->_getTri(i);

        point3d v[] = {mesh->_getVertex(tri[0]),
                       mesh->_getVertex(tri[1]),
                       mesh->_getVertex(tri[2])};

        point3d n = (v[1] - v[0]).cross(v[2] - v[0]);
        n /= std::sqrt(n.dot(n));

        auto m = mesh->getTriNorm(i);
        ASSERT_DOUBLE_EQ(n[0], m[0]);
        ASSERT_DOUBLE_EQ(n[1], m[1]);
        ASSERT_DOUBLE_EQ(n[2], m[2]);
    }
}

TEST_F(TetmeshTest, intersectMontecarlo_all3) {
    double v[] = {0.5, 0.5, 0., 0.5, 0.5, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment all 3 tets should intersect
    ASSERT_EQ(res.size(), 3ul);

    // Intersects tet 0,1,2; exact percentage depends on number of samples
    std::array<double, 3> expected{0.21132, 0.21132, 0.5773};

    // exact percentage depends on number of samples
    double error = 1.0 / nsample;
    for (const auto& r: res)
        ASSERT_LT(std::abs(r.second - expected[r.first.get()]), error);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);

    // Batch mode
    auto res2 = mesh->intersect(v, 2, nsample);
    for (const auto& t: res2)
        for (const auto& r: t)
            std::cout << r.first << ": " << r.second << std::endl;
    ASSERT_EQ(res2.size(), 1ul);
    ASSERT_EQ(res2[0].size(), 3ul);
    sum = 0.0;
    for (const auto& r: res2[0]) {
        ASSERT_LT(std::abs(r.second - expected[r.first]), error);
        sum += r.second;
    }
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectMontecarlo_parallel2face_normal2edge_2) {
    double v[] = {0.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment only 2 tets should intersect
    ASSERT_EQ(res.size(), 2ul);

    // Intersects tet 0,1
    std::array<double, 3> expected{0.5, 0.5, 0.0};

    // exact percentage depends on number of samples
    double error = 1.0 / nsample;
    for (const auto& r: res)
        ASSERT_LT(std::abs(r.second - expected[r.first.get()]), error);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectMontecarlo_parallel2face_normal2edge_3) {
    double v[] = {0.0, 0.0, 0.0, 0.5, 0.866025403784, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment only 2 tets should intersect
    ASSERT_EQ(res.size(), 2ul);

    // Intersects tet 0,2
    std::array<double, 3> expected{0.5, 0.0, 0.5};

    // exact percentage depends on number of samples
    double error = 1.0 / nsample;
    for (const auto& r: res)
        ASSERT_LT(std::abs(r.second - expected[r.first.get()]), error);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectMontecarlo_parallel2face_parallel2edge) {
    double v[] = {0.5, 0.866025403784, 0.0, 0.0, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment intersect only the first tet
    ASSERT_EQ(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectMontecarlo_parallel2face_withinFace) {
    double v[] = {0.2, 0.1, 0.0, 0.5, 0.6, 0.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment can intersect 1 or 2 tets
    ASSERT_GE(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectMontecarlo_parallel2internalFace) {
    double v[] = {0.5, 0.866025403784, 0.0, 0.5, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    constexpr size_t nsample = 100;
    auto res = mesh->intersectMontecarlo<tetrahedron_global_id>(p0, p1, std::nullopt, nsample);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment can intersect 1 or 2 tets
    ASSERT_GE(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_all3_touchingFaces) {
    double v[] = {0.5, 0.5, 0.0, 0.5, 0.5, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment all 3 tets should intersect
    ASSERT_EQ(res.size(), 3ul);

    // Intersects tet 0,1,2
    std::array<double, 3> expected{0.21132, 0.21132, 0.5773};

    for (const auto& r: res)
        ASSERT_GE(r.second, expected[r.first.get()]);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);

    // Batch mode
    auto res2 = mesh->intersect(v, 2);
    for (const auto& t: res2)
        for (const auto& r: t)
            std::cout << r.first << ": " << r.second << std::endl;
    ASSERT_EQ(res2.size(), 1ul);
    ASSERT_EQ(res2[0].size(), 3ul);
    sum = 0.0;
    for (const auto& r: res2[0]) {
        ASSERT_GE(r.second, expected[r.first]);
        sum += r.second;
    }
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_all3_inside) {
    double small = 1e-8;
    double v[] = {0.5, 0.5, small, 0.5, 0.5, 1.0 - small};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment all 3 tets should intersect
    ASSERT_EQ(res.size(), 3ul);

    // Intersects tet 0,1,2
    std::array<double, 3> expected{0.21132, 0.21132, 0.5773};

    for (const auto& r: res)
        ASSERT_GE(r.second, expected[r.first.get()]);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_parallel2face_normal2edge_2) {
    double v[] = {0.0, 0.0, 0.0, 1.0, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment only 2 tets should intersect
    ASSERT_EQ(res.size(), 2ul);

    ASSERT_DOUBLE_EQ(res[0].second, 0.50);
    ASSERT_DOUBLE_EQ(res[1].second, 0.50);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_parallel2face_normal2edge_3) {
    double v[] = {0.0, 0.0, 0.0, 0.5, 0.866025403784, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment only 2 tets should intersect
    ASSERT_EQ(res.size(), 2ul);

    // Intersects tet 0,2
    std::array<double, 3> expected{0.5, 0.0, 0.5};

    for (const auto& r: res)
        ASSERT_DOUBLE_EQ(r.second, expected[r.first.get()]);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_parallel2face_parallel2edge) {
    double v[] = {0.5, 0.866025403784, 0.0, 0.0, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment intersect only the first tet
    ASSERT_EQ(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_parallel2face_withinFace) {
    double v[] = {0.2, 0.1, 0.0, 0.5, 0.6, 0.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment can intersect 1 or 2 tets
    ASSERT_GE(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_parallel2internalFace) {
    double v[] = {0.5, 0.866025403784, 0.0, 0.5, 0.0, 1.0};
    point3d p0(v[0], v[1], v[2]), p1(v[3], v[4], v[5]);

    auto res = mesh->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment can intersect 1 or 2 tets
    ASSERT_GE(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}

TEST_F(TetmeshTest, intersectDeterministic_notFirstNeighbor) {
    // Data from a bug in the field
    std::vector<double> mesh_verteces{
        0.00012257076822803781, 0.0006780064542649207,  0.00042064630907423571,
        0.00018535568938843539, 0.00073659118484427756, 0.00067122005132629644,
        0.00029303616089884128, 0.00045726275951556069, 0.0005383382210917974,
        0.00036866902672842808, 0.00070347379902898656, 0.00045565876460582208,
        9.9156869667801759e-05, 0.00050474632933539439, 0.00074120510599487942,
        -0.0001278801110164713, 0.00044968006878609828, 0.00035668316388897281,
        0.0001655399527923312,  0.00026864657156772679, 0.00025716312239622502,
        0.0001139372139678326,  0.00055231248254049748, 0.0001332087041703196,
        0.00032849825384428558, 0.00051554680707654522, 0.00028303284178744602};
    std::vector<index_t> mesh_tets{0, 1, 2, 3, 0, 2, 1, 4, 0, 5, 2, 4, 0, 5, 6, 2,
                                   6, 0, 7, 5, 0, 6, 7, 8, 0, 2, 6, 8, 0, 2, 8, 3};

    std::unique_ptr<Tetmesh> mesh2(new Tetmesh(mesh_verteces, mesh_tets));


    point3d p0(0.00028457536706320819, 0.00065831720180501655, 0.00045798374339401439),
        p1(0.00028442636498801285, 0.00065874121303548521, 0.00045768874522506909);

    auto res = mesh2->intersectDeterministic<tetrahedron_global_id>(p0, p1);
    for (const auto& r: res)
        std::cout << r.first << ": " << r.second << std::endl;

    // Given the mesh and the segment can intersect 1 or 2 tets
    ASSERT_GE(res.size(), 1ul);

    // sum should be 100%
    double sum{0.0};
    for (const auto& r: res)
        sum += r.second;
    ASSERT_DOUBLE_EQ(sum, 1.0);
}