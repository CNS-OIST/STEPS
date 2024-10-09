

#include <catch2/catch_test_macros.hpp>
#include <catch2/matchers/catch_matchers_floating_point.hpp>

#include "test_common.hpp"

#include "steps/geom/dist/distmesh.hpp"
#include "steps/util/debug.hpp"
#include "steps/util/finish.hpp"
#include "steps/util/init.hpp"

namespace osh = Omega_h;

namespace steps::dist {

using nextIntersectionRes = std::pair<DistMesh::intersectionInfo, DistMesh::intersectionID>;
using point3d = steps::math::point3d;
constexpr auto tol = 4.0 * std::numeric_limits<point3d::value_type>::epsilon();
constexpr auto eps = tol / 100.0;

template <class T>
std::vector<int> conv(const T& v) {
    std::vector<int> ans;
    ans.reserve(v.size());
    for (const auto& i: v) {
        ans.push_back(i.get());
    }
    return ans;
}

bool compare_info(const nextIntersectionRes& a, const nextIntersectionRes& b) {
    return a.second == b.second && a.first.almostEqual(b.first, tol);
}

auto lib = Omega_h::Library();

TEST_CASE("distmesh_basic_access") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    std::vector<int> expected;

    const mesh::tetrahedron_local_id_t tet0{0};
    expected = {2, 1};
    REQUIRE(conv(mesh.getTetTetNeighb(tet0, false)) == expected);
    expected = {3, 0, 7, 1};
    REQUIRE(conv(mesh.getTetTriNeighb(tet0)) == expected);
    expected = {1, 8, 2, 0, 5, 6};
    REQUIRE(conv(mesh.getTetBarNeighb(tet0)) == expected);
    expected = {0, 2, 3, 1};
    REQUIRE(conv(mesh.getTetVertNeighb(tet0)) == expected);

    const mesh::triangle_local_id_t tri0{0};
    expected = {0, 1};
    REQUIRE(conv(mesh.getTriTetNeighb(tri0, false)) == expected);
    expected = {3, 4, 5, 1, 2, 7, 8};
    REQUIRE(conv(mesh.getTriTriNeighbs(tri0, false)) == expected);
    expected = {1, 0, 5};
    REQUIRE(conv(mesh.getTriBarNeighb(tri0)) == expected);
    expected = {2, 0, 1};
    REQUIRE(conv(mesh.getTriVertNeighb(tri0)) == expected);

    const mesh::bar_local_id_t bar0{0};
    expected = {0, 1};
    REQUIRE(conv(mesh.getBarTetNeighb(bar0, false)) == expected);
    expected = {0, 1, 2};
    REQUIRE(conv(mesh.getBarTriNeighb(bar0)) == expected);
    expected = {0, 1};
    REQUIRE(conv(mesh.getBarVertNeighb(bar0)) == expected);

    const mesh::vertex_local_id_t vert0{0};
    expected = {0, 1, 2};
    REQUIRE(conv(mesh.getVertTetNeighb(vert0, false)) == expected);
    expected = {0, 1, 2, 3, 4, 5, 6};
    REQUIRE(conv(mesh.getVertTriNeighb(vert0)) == expected);
    expected = {0, 1, 2, 3, 4};
    REQUIRE(conv(mesh.getVertBarNeighb(vert0)) == expected);
}

TEST_CASE("distmesh_findIntersection") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    DistMesh::intersectionID expected;
    expected = mesh::tetrahedron_local_id_t{2};
    REQUIRE((mesh.findIntersection({-0.1, 0.1, 0.1}).intersection_ == expected));
    expected = mesh::tetrahedron_local_id_t{};
    REQUIRE((mesh.findIntersection({-1, 1, 1}).intersection_ == expected));
    expected = mesh::triangle_local_id_t{5};
    REQUIRE((mesh.findIntersection({-0.1, 0.1, 0}).intersection_ == expected));
    expected = mesh::bar_local_id_t{4};
    REQUIRE((mesh.findIntersection({-0.1, 0, 0}).intersection_ == expected));
    expected = mesh::vertex_local_id_t{3};
    REQUIRE((mesh.findIntersection({0, 0, 1}).intersection_ == expected));
    expected = mesh::vertex_local_id_t{4};
    REQUIRE((mesh.findIntersection({0, 0, -1 + eps}).intersection_ == expected));
    // out of the mesh by eps. Still finds it
    expected = mesh::vertex_local_id_t{4};
    REQUIRE((mesh.findIntersection({0, 0, -1 - eps}).intersection_ == expected));
}

/// check intersectionInfo basic functions
TEST_CASE("distmesh_intersectionInfo") {
    DistMesh::intersectionInfo inIn0{{0, 0, 0}, mesh::tetrahedron_local_id_t{0}};
    DistMesh::intersectionInfo inIn1{{0, 1, 0}, mesh::tetrahedron_local_id_t{0}};
    DistMesh::intersectionInfo inIn2{{0, 0, 0}, mesh::tetrahedron_local_id_t{1}};
    DistMesh::intersectionInfo inIn3{{0, 0, 0}, mesh::triangle_local_id_t{0}};
    DistMesh::intersectionInfo inIn4{{0, eps, 0}, mesh::tetrahedron_local_id_t{0}};

    REQUIRE(inIn0.almostEqual(inIn0, tol));
    REQUIRE(!inIn0.almostEqual(inIn1, tol));
    REQUIRE(!inIn0.almostEqual(inIn2, tol));
    REQUIRE(!inIn0.almostEqual(inIn3, tol));
    REQUIRE(inIn0.almostEqual(inIn4, tol));
}

TEST_CASE("distmesh_findNextIntersection_fromTet") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    // p_beg in tet
    const auto p_beg_info = mesh.findIntersection({0.1, 0.1, -0.1});
    // p_end in tet
    auto p_end_info = mesh.findIntersection({0.1, 0.2, -0.1});
    nextIntersectionRes expected{{{0.1, 0.2, -0.1}, mesh::tetrahedron_local_id_t{1}},
                                 mesh::tetrahedron_local_id_t{1}};
    auto res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end in another tet
    p_end_info = mesh.findIntersection({0.1, 0.1, 0.1});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on vertex
    p_end_info = mesh.findIntersection({0, 0, -1});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh. We should receive the vertex instead
    p_end_info = mesh.findIntersection({0, 0, -1 - eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely inside the mesh. return the vertex or something very close to it
    p_end_info = mesh.findIntersection({0, 0, -1 + eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh and along the direction of the vertex. Return the vertex again
    p_end_info = mesh.findIntersection({-0.1, -0.1, -2 + 0.1});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on bar
    p_end_info = mesh.findIntersection({0, 0, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet. same intersection on bar
    p_end_info = mesh.findIntersection({eps, eps, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet. same intersection on bar
    p_end_info = mesh.findIntersection({-eps, -eps, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of tet. same intersection on bar
    p_end_info = mesh.findIntersection({-0.1, -0.1, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on tri
    p_end_info = mesh.findIntersection({0.1, 0, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet. same intersection on bar
    p_end_info = mesh.findIntersection({0.1, eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet. same intersection on bar
    p_end_info = mesh.findIntersection({0.1, -eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of tet. same intersection on bar
    p_end_info = mesh.findIntersection({0.1, -0.1, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
}

TEST_CASE("distmesh_findNextIntersection_fromTri") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    // p_beg on tri
    const auto p_beg_info = mesh.findIntersection({0.1, 0, -0.1});
    // p_end in tet
    auto p_end_info = mesh.findIntersection({0.1, 0.1, -0.1});
    nextIntersectionRes expected{p_end_info, mesh::tetrahedron_local_id_t{1}};
    auto res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on outer tri as well
    p_end_info = mesh.findIntersection({0.2, 0, -0.1});
    expected = nextIntersectionRes{mesh.findIntersection({0.2, 0, -0.1}),
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet. same result
    p_end_info = mesh.findIntersection({0.2, eps, -0.1});
    expected = nextIntersectionRes{mesh.findIntersection({0.2, 0, -0.1}),
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet. same result
    p_end_info = mesh.findIntersection({0.2, -eps, -0.1});
    expected = nextIntersectionRes{mesh.findIntersection({0.2, 0, -0.1}),
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on vert
    p_end_info = mesh.findIntersection({0, 0, -1});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet. same intersection
    p_end_info = mesh.findIntersection({eps, 0, -1 + eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({eps, eps, -1 + eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet. same intersection
    p_end_info = mesh.findIntersection({0, 0, -1 - eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of tet. vert on segment
    p_end_info = mesh.findIntersection({-0.1, 0, -1.9});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on shared edge
    p_end_info = mesh.findIntersection({0.1, 0, 0});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tri
    p_end_info = mesh.findIntersection({0.1, 0, -eps});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tri but in mesh
    p_end_info = mesh.findIntersection({0.1, 0, eps});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on another tri
    p_end_info = mesh.findIntersection({0.1, 0, 0.1});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on boundary edge
    p_end_info = mesh.findIntersection({0, 0, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tri
    p_end_info = mesh.findIntersection({eps, 0, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh
    p_end_info = mesh.findIntersection({-eps, 0, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh
    p_end_info = mesh.findIntersection({-0.1, 0, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on another boundary tri
    p_end_info = mesh.findIntersection({0, 0.1, -0.1});
    expected = nextIntersectionRes{{{0, 0.1, -0.1}, mesh::triangle_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({eps, 0.1, -0.1});
    expected = nextIntersectionRes{{{0, 0.1, -0.1}, mesh::triangle_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet
    p_end_info = mesh.findIntersection({-eps, 0.1, -0.1});
    expected = nextIntersectionRes{{{0, 0.1, -0.1}, mesh::triangle_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of tet
    p_end_info = mesh.findIntersection({-0.1, 0.2, -0.1});
    expected = nextIntersectionRes{{{0, 0.1, -0.1}, mesh::triangle_local_id_t{4}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on shared tri
    p_end_info = mesh.findIntersection({0.1, 0.1, 0});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0.1, 0.1, -eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in another tet
    p_end_info = mesh.findIntersection({0.1, 0.1, eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end in another tet
    p_end_info = mesh.findIntersection({0.1, 0.2, 0.1});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on opposite vert
    p_end_info = mesh.findIntersection({0, 1, 0});
    expected = nextIntersectionRes{{{0, 1, 0}, mesh::vertex_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0, 1 - eps, 0});
    expected = nextIntersectionRes{{{0, 1, 0}, mesh::vertex_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh
    p_end_info = mesh.findIntersection({0, 1 - eps, 0});
    expected = nextIntersectionRes{{{0, 1, 0}, mesh::vertex_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in another tet
    p_end_info = mesh.findIntersection({0, 1 - eps, eps / 2.0});
    expected = nextIntersectionRes{{{0, 1, 0}, mesh::vertex_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh
    p_end_info = mesh.findIntersection({-0.1, 2, 0.1});
    expected = nextIntersectionRes{{{0, 1, 0}, mesh::vertex_local_id_t{2}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on bar, not tri bar
    p_end_info = mesh.findIntersection({0, 0.1, 0});
    expected = nextIntersectionRes{{{0, 0.1, 0}, mesh::bar_local_id_t{1}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({eps, 0.1, -eps});
    expected = nextIntersectionRes{{{0, 0.1, 0}, mesh::bar_local_id_t{1}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in another tets
    p_end_info = mesh.findIntersection({-eps, 0.1, eps});
    expected = nextIntersectionRes{{{0, 0.1, 0}, mesh::bar_local_id_t{1}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh
    p_end_info = mesh.findIntersection({-eps, 0.1, -eps});
    expected = nextIntersectionRes{{{0, 0.1, 0}, mesh::bar_local_id_t{1}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end in another tet
    p_end_info = mesh.findIntersection({-0.1, 0.2, 0.1});
    expected = nextIntersectionRes{{{0, 0.1, 0}, mesh::bar_local_id_t{1}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
}

TEST_CASE("distmesh_findNextIntersection_fromBar") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    // p_beg along bar
    const auto p_beg_info = mesh.findIntersection({0, 0, -0.1});
    // p_end in tet
    auto p_end_info = mesh.findIntersection({0.1, 0.1, -0.1});
    nextIntersectionRes expected{mesh.findIntersection({0.1, 0.1, -0.1}),
                                 mesh::tetrahedron_local_id_t{1}};
    auto res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on adjacent tri
    p_end_info = mesh.findIntersection({0.1, 0, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0.1, eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh
    p_end_info = mesh.findIntersection({0.1, -eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on bar
    p_end_info = mesh.findIntersection({0, 0, -0.2}),
    expected = nextIntersectionRes{{{0, 0, -0.2}, mesh::bar_local_id_t{3}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on connected vertex
    p_end_info = mesh.findIntersection({0, 0, -1});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end close to vertex
    p_end_info = mesh.findIntersection({0, 0, -1 + eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0, 0, -1 - eps});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of bar. vert on segment
    p_end_info = mesh.findIntersection({0, 0, -2});
    expected = nextIntersectionRes{{{0, 0, -1}, mesh::vertex_local_id_t{4}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on opposite bar
    p_end_info = mesh.findIntersection({0.5, 0.5, 0});
    expected = nextIntersectionRes{{{0.5, 0.5, 0}, mesh::bar_local_id_t{5}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely outside
    p_end_info = mesh.findIntersection({0.5, 0.5 + eps, 0});
    expected = nextIntersectionRes{{{0.5, 0.5, 0}, mesh::bar_local_id_t{5}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely inside
    p_end_info = mesh.findIntersection({0.5 - eps, 0.5 - eps, -eps / 2});
    expected = nextIntersectionRes{{{0.5, 0.5, 0}, mesh::bar_local_id_t{5}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end outside mesh
    p_end_info = mesh.findIntersection({1, 1, 0.1});
    expected = nextIntersectionRes{{{0.5, 0.5, 0}, mesh::bar_local_id_t{5}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on vert on connected tri
    p_end_info = mesh.findIntersection({1, 0, 0});
    expected = nextIntersectionRes{{{1, 0, 0}, mesh::vertex_local_id_t{1}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on vert on connected tri, barely out
    p_end_info = mesh.findIntersection({1 + eps, 0, 0});
    expected = nextIntersectionRes{{{1, 0, 0}, mesh::vertex_local_id_t{1}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on vert on connected tri, barely in
    p_end_info = mesh.findIntersection({1 - eps, 0, 0});
    expected = nextIntersectionRes{{{1, 0, 0}, mesh::vertex_local_id_t{1}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh. vert on connected tri is the intersection
    p_end_info = mesh.findIntersection({2, 0, 0.1});
    expected = nextIntersectionRes{{{1, 0, 0}, mesh::vertex_local_id_t{1}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on vert-connected try (to the bar)
    p_end_info = mesh.findIntersection({0.1, 0.1, 0});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0.1, 0.1, -eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0.1, 0.1, eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0.2, 0.2, 0.1});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
}

TEST_CASE("distmesh_findNextIntersection_fromVert") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    // p_beg on vert
    const auto p_beg_info = mesh.findIntersection({0, 0, -1});

    // p_end on connected bar
    auto p_end_info = mesh.findIntersection({0, 0, -0.1});
    nextIntersectionRes expected{{{0, 0, -0.1}, mesh::bar_local_id_t{3}}, mesh::bar_local_id_t{3}};
    auto res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in mesh
    p_end_info = mesh.findIntersection({eps, eps, -0.1 + eps});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of mesh
    p_end_info = mesh.findIntersection({-eps, -eps, -0.1});
    expected = nextIntersectionRes{{{0, 0, -0.1}, mesh::bar_local_id_t{3}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on connected vert
    p_end_info = mesh.findIntersection({0, 0, 0});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0, 0, -eps});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely out of tet
    p_end_info = mesh.findIntersection({0, 0, eps});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end in another tet
    p_end_info = mesh.findIntersection({0, 0, 0.1});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh and segment along bar
    p_end_info = mesh.findIntersection({0, 0, 2});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end in tet but initial part of segment out of tet
    p_end_info = mesh.findIntersection({-0.1, 0.1, 0.1});
    expected = nextIntersectionRes{{p_beg_info.point_, mesh::tetrahedron_local_id_t{}},
                                   mesh::tetrahedron_local_id_t{}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    REQUIRE(mesh.findIntersection({-0.1, 0.1, 0.1})
                .almostEqual(DistMesh::intersectionInfo{{-0.1, 0.1, 0.1},
                                                        mesh::tetrahedron_local_id_t{2}},
                             tol));

    // p_end on connected tri
    p_end_info = mesh.findIntersection({0.1, 0, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0.1, eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh
    p_end_info = mesh.findIntersection({0.1, -eps, -0.1});
    expected = nextIntersectionRes{{{0.1, 0, -0.1}, mesh::triangle_local_id_t{2}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on opposite tri
    p_end_info = mesh.findIntersection({0.1, 0.1, 0});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on opposite tri, barely out
    p_end_info = mesh.findIntersection({0.1, 0.1, eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end on opposite tri, barely in
    p_end_info = mesh.findIntersection({0.1, 0.1, -eps});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh, intersection on tri
    p_end_info = mesh.findIntersection({0.2, 0.2, 1});
    expected = nextIntersectionRes{{{0.1, 0.1, 0}, mesh::triangle_local_id_t{0}},
                                   mesh::tetrahedron_local_id_t{1}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end on opposite bar
    p_end_info = mesh.findIntersection({0.1, 0, 0});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in tet
    p_end_info = mesh.findIntersection({0.1, eps, -eps});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end barely in other tet
    p_end_info = mesh.findIntersection({0.1, eps, eps});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh. segment - bar intersection
    p_end_info = mesh.findIntersection({0.2, 0, 1});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    // p_end out of mesh
    p_end_info = mesh.findIntersection({-eps, -eps, 1});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));

    // p_end barely out of mesh in various ways
    p_end_info = mesh.findIntersection({-eps, -eps, -eps});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({-eps, 0, -eps});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0, -eps, -eps});
    expected = nextIntersectionRes{{{0, 0, 0}, mesh::vertex_local_id_t{0}},
                                   mesh::bar_local_id_t{3}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
    p_end_info = mesh.findIntersection({0.1, -eps, 0});
    expected = nextIntersectionRes{{{0.1, 0, 0}, mesh::bar_local_id_t{0}},
                                   mesh::triangle_local_id_t{2}};
    res = mesh.findNextIntersection(p_beg_info, p_end_info);
    CAPTURE(p_beg_info, p_end_info, res, expected);
    REQUIRE(compare_info(res, expected));
}

#define REQUIRE_EQUAL_INTERSECTION_LISTS(ans, required)                      \
    CAPTURE(ans, required);                                                  \
    do {                                                                     \
        REQUIRE(ans.size() == required.size());                              \
        for (std::size_t i = 0; i < ans.size(); ++i) {                       \
            REQUIRE(ans[i].size() == required[i].size());                    \
            for (std::size_t j = 0; j < ans[i].size(); ++j) {                \
                REQUIRE(ans[i][j].first == required[i][j].first);            \
                REQUIRE(std::abs(ans[i][j].second - required[i][j].second) < \
                        steps::math::tol_lin);                               \
            }                                                                \
        }                                                                    \
    } while (0)

TEST_CASE("distmesh_intersect") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    const auto compute_ans = [&](const std::vector<double>& pp) {
        auto ans = mesh.intersect(pp.data(), pp.size() / 3);
        for (auto& i: ans) {
            std::sort(i.begin(), i.end(), [](const auto& p0, const auto& p1) {
                return p0.first < p1.first;
            });
        }
        return ans;
    };

    std::vector<double> pp;
    std::vector<DistMesh::intersection_list_t> expected;

    // // same tet
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, -0.2};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    compute_ans(pp);
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // 2 tets
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, 0.1};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 0.5}, {mesh::tetrahedron_local_id_t{1}, 0.5}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // 2 tets, on the connecting tri
    pp = std::vector<double>{0.1, 0.1, 0, 0.1, 0.2, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 0.5}, {mesh::tetrahedron_local_id_t{1}, 0.5}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // on bar
    pp = std::vector<double>{0, 0, -0.1, 0, 0, -0.2};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // on bar, touching a vertex
    pp = std::vector<double>{0, 0, 0, 0, 0, -0.2};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // on bar, between 2 tets
    pp = std::vector<double>{0.1, 0, 0, 0.2, 0, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 0.5}, {mesh::tetrahedron_local_id_t{1}, 0.5}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // concave mesh and line on tris
    pp = std::vector<double>{0.1, 0, -0.5, -0.5, 0, 0.1};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{1}, 1.0 / 6.0},
         {mesh::tetrahedron_local_id_t{2}, 1.0 / 6.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // on bar, among 3 tets
    pp = std::vector<double>{0, 0.1, 0, 0, 0.2, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 1.0 / 3.0},
         {mesh::tetrahedron_local_id_t{1}, 1.0 / 3.0},
         {mesh::tetrahedron_local_id_t{2}, 1.0 / 3.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // on bar, among 3 tets, starting from vert
    pp = std::vector<double>{0, 0, 0, 0, 0.2, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 1.0 / 3.0},
         {mesh::tetrahedron_local_id_t{1}, 1.0 / 3.0},
         {mesh::tetrahedron_local_id_t{2}, 1.0 / 3.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // starting from vert
    pp = std::vector<double>{0, 0, -1.0, 0.1, 0.1, -0.1};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // starting shared vert
    pp = std::vector<double>{0, 0, 0, 0.1, 0.1, -0.1};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // cross 2 tets and a vert in the middle with more shared tets
    pp = std::vector<double>{-0.1, 0, 0.1, 0.1, 0, -0.1};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{1}, 0.5}, {mesh::tetrahedron_local_id_t{2}, 0.5}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // cross 2 tets and a bar in the middle with more shared tets
    pp = std::vector<double>{-0.1, 0.1, 0.1, 0.1, 0.1, -0.1};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{1}, 0.5}, {mesh::tetrahedron_local_id_t{2}, 0.5}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // go out of the mesh. The intersection point is at {1/3, 1/3, 1/3}
    pp = std::vector<double>{0.1, 0.1, 0.1, 0.9, 0.9, 0.9};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 7.0 / 24.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // go out of the mesh on a tri
    pp = std::vector<double>{0.1, 0.1, 0, 0.9, 0.9, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 0.25}, {mesh::tetrahedron_local_id_t{1}, 0.25}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // go out of the mesh on a shared bar starting from a shared vert
    pp = std::vector<double>{0, 0, 0, 0, 2, 0};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 1.0 / 6.0},
         {mesh::tetrahedron_local_id_t{1}, 1.0 / 6.0},
         {mesh::tetrahedron_local_id_t{2}, 1.0 / 6.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // from out of mesh
    pp = std::vector<double>{0.9, 0.9, 0.9, 0.1, 0.1, 0.1};
    expected = std::vector<DistMesh::intersection_list_t>{
        {{mesh::tetrahedron_local_id_t{0}, 7.0 / 24.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // always out of mesh
    pp = std::vector<double>{0.9, 0.9, 0.9, 1, 1, 1};
    expected = std::vector<DistMesh::intersection_list_t>{{}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // always out of mesh but a vert touces
    pp = std::vector<double>{1, 0, 0, 2, 0, 0};
    expected = std::vector<DistMesh::intersection_list_t>{{}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // always out of mesh but a vert touces
    pp = std::vector<double>{-1, 0, -1, 1, 0, -1};
    expected = std::vector<DistMesh::intersection_list_t>{{}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // intersect 2 segments
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, -0.2, 0.1, 0.1, -0.3};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}},
                                                          {{mesh::tetrahedron_local_id_t{1}, 1.0}}};

    // intersect 2 segments in 2 tets
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, 0, 0.1, 0.1, 0.1};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}},
                                                          {{mesh::tetrahedron_local_id_t{0}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);
}

TEST_CASE("distmesh_intersectIndependentSegments") {
    const auto mesh_file = Omega_h::filesystem::path(STEPS_SOURCE_DIR) / "test" / "mesh" /
                           "3_tets.msh";
    steps::dist::DistMesh mesh(lib, mesh_file.string());

    const auto compute_ans = [&](const std::vector<double>& pp) {
        auto ans = mesh.intersectIndependentSegments(pp.data(), pp.size() / 3);
        for (auto& i: ans) {
            std::sort(i.begin(), i.end(), [](const auto& p0, const auto& p1) {
                return p0.first < p1.first;
            });
        }
        return ans;
    };

    std::vector<double> pp;
    std::vector<DistMesh::intersection_list_t> expected;

    // // same tet
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, -0.2};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}}};
    compute_ans(pp);
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);

    // intersect 2 segments in 2 tets
    pp = std::vector<double>{0.1, 0.1, -0.1, 0.1, 0.1, 0, 0.1, 0.1, 0, 0.1, 0.1, 0.1};
    expected = std::vector<DistMesh::intersection_list_t>{{{mesh::tetrahedron_local_id_t{1}, 1.0}},
                                                          {{mesh::tetrahedron_local_id_t{0}, 1.0}}};
    REQUIRE_EQUAL_INTERSECTION_LISTS(compute_ans(pp), expected);
}
}  // namespace steps::dist
