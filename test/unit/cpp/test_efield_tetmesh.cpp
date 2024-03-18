#include <cmath>

#include "solver/efield/tetmesh.hpp"

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

using namespace steps::solver::efield;

TEST_CASE("Efield_setTriCapac") {
    using steps::vertex_id_t;
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
    vertex_id_t tris[12]{vertex_id_t(0),
                         vertex_id_t(1),
                         vertex_id_t(2),
                         vertex_id_t(0),
                         vertex_id_t(1),
                         vertex_id_t(3),
                         vertex_id_t(1),
                         vertex_id_t(2),
                         vertex_id_t(3),
                         vertex_id_t(0),
                         vertex_id_t(2),
                         vertex_id_t(3)};
    vertex_id_t tets[4]{vertex_id_t(0), vertex_id_t(1), vertex_id_t(2), vertex_id_t(3)};

    TetMesh mesh{nverts, verts, ntris, tris, ntets, tets};

    const double triArea{std::sqrt(3.0)};

    mesh.allocateSurface();

    REQUIRE(mesh.getTotalCapacitance() == 0);
    for (auto i: steps::vertex_id_t::range(mesh.countVertices())) {
        REQUIRE(mesh.getVertex(i)->getCapacitance() == 0);
    }

    mesh.applySurfaceCapacitance(1.0);
    REQUIRE_THAT(mesh.getTotalCapacitance(), Catch::Matchers::WithinAbs(4.0 * triArea, abs_err));
    for (auto i: steps::vertex_id_t::range(mesh.countVertices())) {
        REQUIRE_THAT(mesh.getVertex(i)->getCapacitance(),
                     Catch::Matchers::WithinAbs(triArea, abs_err));
    }

    mesh.applyTriCapacitance(steps::triangle_local_id(0), 2.0);
    REQUIRE_THAT(mesh.getTotalCapacitance(), Catch::Matchers::WithinAbs(5.0 * triArea, abs_err));
    for (auto i: steps::vertex_id_t::range(3)) {
        REQUIRE_THAT(mesh.getVertex(i)->getCapacitance(),
                     Catch::Matchers::WithinAbs(4.0 * triArea / 3.0, abs_err));
    }
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(3))->getCapacitance(),
                 Catch::Matchers::WithinAbs(triArea, abs_err));

    mesh.applyTriCapacitance(steps::triangle_local_id(2), 3.0);
    REQUIRE_THAT(mesh.getTotalCapacitance(), Catch::Matchers::WithinAbs(7.0 * triArea, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(0))->getCapacitance(),
                 Catch::Matchers::WithinAbs(4.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(1))->getCapacitance(),
                 Catch::Matchers::WithinAbs(6.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(2))->getCapacitance(),
                 Catch::Matchers::WithinAbs(6.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(3))->getCapacitance(),
                 Catch::Matchers::WithinAbs(5.0 * triArea / 3.0, abs_err));

    mesh.applyTriCapacitance(steps::triangle_local_id(0), 1.0);
    mesh.applyTriCapacitance(steps::triangle_local_id(2), 1.0);
    REQUIRE_THAT(mesh.getTotalCapacitance(), Catch::Matchers::WithinAbs(4.0 * triArea, abs_err));
    for (auto i: steps::vertex_id_t::range(mesh.countVertices())) {
        REQUIRE_THAT(mesh.getVertex(i)->getCapacitance(),
                     Catch::Matchers::WithinAbs(triArea, abs_err));
    }

    for (auto i: steps::triangle_local_id::range(ntris)) {
        mesh.applyTriCapacitance(i, 0.0);
    }
    REQUIRE(mesh.getTotalCapacitance() == 0);
    for (auto i: steps::vertex_id_t::range(mesh.countVertices())) {
        REQUIRE(mesh.getVertex(i)->getCapacitance() == 0);
    }

    for (auto i: steps::triangle_local_id::range(ntris)) {
        mesh.applyTriCapacitance(i, i.get());
    }
    REQUIRE_THAT(mesh.getTotalCapacitance(), Catch::Matchers::WithinAbs(6.0 * triArea, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(0))->getCapacitance(),
                 Catch::Matchers::WithinAbs(4.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(1))->getCapacitance(),
                 Catch::Matchers::WithinAbs(3.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(2))->getCapacitance(),
                 Catch::Matchers::WithinAbs(5.0 * triArea / 3.0, abs_err));
    REQUIRE_THAT(mesh.getVertex(steps::vertex_id_t(3))->getCapacitance(),
                 Catch::Matchers::WithinAbs(6.0 * triArea / 3.0, abs_err));
}
