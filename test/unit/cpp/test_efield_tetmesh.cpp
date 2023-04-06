#include <cmath>

#include "solver/efield/tetmesh.hpp"

#include <gtest/gtest.h>

using namespace steps::solver::efield;

TEST(Efield, setTriCapac) {
    double abs_err{1e-10};

    uint nverts{4}, ntris{4}, ntets{1};
    double verts[12]{1,
                     0,
                     -1.0 / std::sqrt(2.0),
                     -1,
                     0,
                     -1.0 / std::sqrt(2.0),
                     0,
                     1,
                     1.0 / std::sqrt(2.0),
                     0,
                     -1,
                     1.0 / std::sqrt(2.0)};
    steps::vertex_id_t tris[12]{0, 1, 2, 0, 1, 3, 1, 2, 3, 0, 2, 3};
    steps::vertex_id_t tets[4]{0, 1, 2, 3};

    TetMesh mesh{nverts, verts, ntris, tris, ntets, tets};

    const double triArea{std::sqrt(3.0)};

    mesh.allocateSurface();

    ASSERT_EQ(mesh.getTotalCapacitance(), 0);
    for (uint i = 0; i < mesh.countVertices(); ++i) {
        ASSERT_EQ(mesh.getVertex(i)->getCapacitance(), 0);
    }

    mesh.applySurfaceCapacitance(1.0);
    ASSERT_NEAR(mesh.getTotalCapacitance(), 4.0 * triArea, abs_err);
    for (uint i = 0; i < mesh.countVertices(); ++i) {
        ASSERT_NEAR(mesh.getVertex(i)->getCapacitance(), triArea, abs_err);
    }

    mesh.applyTriCapacitance(0, 2.0);
    ASSERT_NEAR(mesh.getTotalCapacitance(), 5.0 * triArea, abs_err);
    for (uint i = 0; i < 3; ++i) {
        ASSERT_NEAR(mesh.getVertex(i)->getCapacitance(), 4.0 * triArea / 3.0, abs_err);
    }
    ASSERT_NEAR(mesh.getVertex(3)->getCapacitance(), triArea, abs_err);

    mesh.applyTriCapacitance(2, 3.0);
    ASSERT_NEAR(mesh.getTotalCapacitance(), 7.0 * triArea, abs_err);
    ASSERT_NEAR(mesh.getVertex(0)->getCapacitance(), 4.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(1)->getCapacitance(), 6.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(2)->getCapacitance(), 6.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(3)->getCapacitance(), 5.0 * triArea / 3.0, abs_err);

    mesh.applyTriCapacitance(0, 1.0);
    mesh.applyTriCapacitance(2, 1.0);
    ASSERT_NEAR(mesh.getTotalCapacitance(), 4.0 * triArea, abs_err);
    for (uint i = 0; i < mesh.countVertices(); ++i) {
        ASSERT_NEAR(mesh.getVertex(i)->getCapacitance(), triArea, abs_err);
    }

    for (uint i = 0; i < ntris; ++i) {
        mesh.applyTriCapacitance(i, 0.0);
    }
    ASSERT_EQ(mesh.getTotalCapacitance(), 0);
    for (uint i = 0; i < mesh.countVertices(); ++i) {
        ASSERT_EQ(mesh.getVertex(i)->getCapacitance(), 0);
    }

    for (uint i = 0; i < ntris; ++i) {
        mesh.applyTriCapacitance(i, i);
    }
    ASSERT_NEAR(mesh.getTotalCapacitance(), 6.0 * triArea, abs_err);
    ASSERT_NEAR(mesh.getVertex(0)->getCapacitance(), 4.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(1)->getCapacitance(), 3.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(2)->getCapacitance(), 5.0 * triArea / 3.0, abs_err);
    ASSERT_NEAR(mesh.getVertex(3)->getCapacitance(), 6.0 * triArea / 3.0, abs_err);
}
