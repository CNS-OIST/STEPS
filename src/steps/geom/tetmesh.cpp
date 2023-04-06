/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

#include "tetmesh.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <numeric>
#include <sstream>
#include <string>
#include <unordered_map>

#include <easylogging++.h>

#include "math/smallsort.hpp"
#include "math/tetrahedron.hpp"
#include "math/triangle.hpp"
#include "util/checkid.hpp"
#include "util/collections.hpp"
#include "util/error.hpp"
#include "util/strong_id.hpp"

using steps::util::as_vector;
using steps::util::checkID;
using steps::util::deref_strongid;
using steps::util::deref_strong_id_func;
using steps::util::make_unique_indexer;

////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace tetmesh {

const Tetmesh::tri_tets Tetmesh::UNKNOWN_TRI_NEIGHBORS{
    {std::nullopt, std::nullopt}};
const Tetmesh::tet_tets Tetmesh::UNKNOWN_TET_NEIGHBORS{
    {std::nullopt, std::nullopt, std::nullopt, std::nullopt}};

Tetmesh::Tetmesh(std::vector<double> const &verts,
                 std::vector<index_t> const &tets,
                 std::vector<index_t> const &tris) {
  using std::to_string;
  using steps::math::small_sort;

  srand(time(nullptr));

  // check the vectors are of the expected size

  ArgErrLogIf(verts.size() % 3 || tris.size() % 3 || tets.size() % 4,
              "Tables supplied to Tet mesh initialiser function are not of the "
              "expected dimensions");

  pVertsN = verts.size() / 3;

  ArgErrLogIf(pVertsN == 0u, "Empty vertex list");

  pTetsN = tets.size() / 4;

  ArgErrLogIf(!pTetsN, "Empty tets list");

  // copy vertices and update bounding box
  pVerts.resize(pVertsN);
  for (size_t i = 0, j = 0; i < pVertsN; ++i, j += 3) {
    point3d vert = {verts[j], verts[j + 1], verts[j + 2]};
    pVerts[i] = vert;
    pBBox.insert(vert);
  }

  // copy tets;
  pTets.resize(pTetsN);
  for (size_t i = 0, j = 0; i < pTetsN; ++i, j += 4) {
    pTets[i] = tet_verts{{tets[j], tets[j + 1], tets[j + 2], tets[j + 3]}};
  }

  // Add user-supplied tris and faces for each tet to pTris and set
  // tet->tri adjacency
  {
    pTet_tri_neighbours.resize(pTetsN);
    auto tri_indices =
        make_unique_indexer<tri_verts>(std::back_inserter(pTris));

    // first add user-supplied tris:
    for (size_t i = 0; i < tris.size(); i += 3) {
      tri_indices.insert(
          small_sort<3>(tri_verts{{tris[i], tris[i + 1], tris[i + 2]}}));
    }

    for (uint i = 0; i < pTetsN; ++i) {
      const tet_verts &tet = pTets[i];
      tri_verts tris_sorted[4] = {
          small_sort<3>(tri_verts{{tet[0], tet[1], tet[2]}}),
          small_sort<3>(tri_verts{{tet[0], tet[1], tet[3]}}),
          small_sort<3>(tri_verts{{tet[0], tet[2], tet[3]}}),
          small_sort<3>(tri_verts{{tet[1], tet[2], tet[3]}})};

      for (int j = 0; j < 4; ++j) {
        pTet_tri_neighbours[i][j] = tri_indices[tris_sorted[j]];
      }
    }
  }
  pTris.shrink_to_fit();
  pTrisN = pTris.size();

  // For each tet, compute volume and barycenter, and update
  // tri->tet adjacency and tet->tet adjacency information.

  pTet_vols.resize(pTetsN);
  pTet_barycenters.resize(pTetsN);

  pTet_tet_neighbours.assign(pTetsN, UNKNOWN_TET_NEIGHBORS);
  pTri_tet_neighbours.assign(pTrisN, UNKNOWN_TRI_NEIGHBORS);

  for (auto i = 0u; i < pTetsN; ++i) {
    auto tet = pTets[i];
    point3d v[4] = {pVerts[tet[0].get()], pVerts[tet[1].get()],
                    pVerts[tet[2].get()], pVerts[tet[3].get()]};

    pTet_vols[i] = steps::math::tet_vol(v[0], v[1], v[2], v[3]);
    pTet_barycenters[i] = steps::math::tet_barycenter(v[0], v[1], v[2], v[3]);

    ArgErrLogIf(pTet_vols[i] <= 0, "degenerate tetrahedron " + to_string(i));

    for (auto face = 0u; face < 4; ++face) {
      auto tri = pTet_tri_neighbours[i][face];
      auto &tri_tets_neigh = pTri_tet_neighbours[tri.get()];

      if (tri_tets_neigh[0].unknown()) {
        tri_tets_neigh[0] = i;
      } else if (tri_tets_neigh[1].unknown()) {
        tri_tets_neigh[1] = i;

        // tri_tets_neigh[0] and [1] are neighbours.
        auto other_tet = tri_tets_neigh[0];
        pTet_tet_neighbours[i][face] = other_tet;

        bool found_other_face = false;
        for (int other_face = 0; other_face < 4; ++other_face) {
          if (pTet_tri_neighbours[other_tet.get()][other_face] == tri) {
              ProgErrLogIf(pTet_tet_neighbours[other_tet.get()][other_face].valid(),
                           "inconsistent tet<->tet association");

              pTet_tet_neighbours[other_tet.get()][other_face] = i;
              found_other_face = true;
              break;
          }
        }

        ProgErrLogIf(!found_other_face, "inconsistent tet<->tet association");

      } else {
        ProgErrLog("inconsisent tri<->tet association");
      }
    }
  }

  // for each tri, compute area, barycenter, norm; allocate vectors for patches
  // and diff. boundaries.

  pTri_areas.resize(pTrisN);
  pTri_barycs.resize(pTrisN);
  pTri_norms.resize(pTrisN);

  pTri_diffboundaries.assign(pTrisN, nullptr);
  pTri_patches.assign(pTrisN, nullptr);

  for (auto i = 0u; i < pTrisN; ++i) {
    auto tri = pTris[i];
    point3d v[3] = {pVerts[tri[0].get()], pVerts[tri[1].get()],
                    pVerts[tri[2].get()]};

    pTri_areas[i] = steps::math::tri_area(v[0], v[1], v[2]);
    pTri_norms[i] = steps::math::tri_normal(v[0], v[1], v[2]);
    pTri_barycs[i] = steps::math::tri_barycenter(v[0], v[1], v[2]);

    ArgErrLogIf(pTri_areas[i] <= 0, "degenerate triangle " + to_string(i));
  }

  // initialise tet compartment association.
  pTet_comps.assign(pTetsN, nullptr);

  // create tri->bar adjacency
  buildBarData();
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::buildBarData() {
  using steps::math::small_sort;

  pTri_bars.resize(pTrisN);

  pBars.clear();
  auto bar_indices = make_unique_indexer<bar_verts>(std::back_inserter(pBars));

  for (uint i = 0; i < pTrisN; ++i) {
    const tri_verts &tri = pTris[i];
    bar_verts bars[3] = {
        small_sort<2>(bar_verts{{tri[0], tri[1]}}),
        small_sort<2>(bar_verts{{tri[0], tri[2]}}),
        small_sort<2>(bar_verts{{tri[1], tri[2]}}),
    };

    for (int j = 0; j < 3; ++j) {
      pTri_bars[i][j] = bar_indices[bars[j]];
    }
  }
  pBars.shrink_to_fit();
  pBarsN = pBars.size();

  // Include surface diffusion boundary stuff here
  pBar_sdiffboundaries.assign(pBarsN, nullptr);

  pBar_tri_neighbours.resize(pBarsN);
  pBar_tri_neighbours.assign(pBarsN, bar_tris{{std::nullopt, std::nullopt}});
}

////////////////////////////////////////////////////////////////////////////////

Tetmesh::Tetmesh(std::vector<double> const &verts,
                 std::vector<vertex_id_t> const &tris,
                 std::vector<double> const &tri_areas,
                 std::vector<double> const &tri_norms,
                 std::vector<tetrahedron_id_t> const &tri_tet_neighbs,
                 std::vector<vertex_id_t> const &tets,
                 std::vector<double> const &tet_vols,
                 std::vector<double> const &tet_barycs,
                 std::vector<triangle_id_t> const &tet_tri_neighbs,
                 std::vector<tetrahedron_id_t> const &tet_tet_neighbs)
    : pVertsN(), pBarsN(0), pTrisN(0), pTetsN(0), pMembs(), pDiffBoundaries() {
  // check the vectors are of the expected size

  ArgErrLogIf(verts.size() % 3 || tris.size() % 3 || tets.size() % 4,
              "Tables supplied to Tet mesh initialiser function are not of the "
              "expected dimensions");

  pVertsN = verts.size() / 3;

  ArgErrLogIf(pVertsN == 0u, "Empty vertex list");

  pTetsN = tets.size() / 4;

  ArgErrLogIf(!pTetsN, "Empty tets list");

  pTrisN = tris.size() / 3;

  ArgErrLogIf(!pTetsN, "Empty tris list");

  ArgErrLogIf(tri_areas.size() != pTrisN, "Inconsistent tri_areas size");

  ArgErrLogIf(tri_norms.size() != pTrisN * 3, "Inconsistent tri_norms size");

  ArgErrLogIf(tri_tet_neighbs.size() != pTrisN * 2,
              "Inconsistent tri_tet_neighbrs size");

  ArgErrLogIf(tet_vols.size() != pTetsN, "Inconsistent tet_vols size");

  ArgErrLogIf(tet_barycs.size() != pTetsN * 3, "Inconsistent tet_barycs size");

  ArgErrLogIf(tet_tri_neighbs.size() != pTetsN * 4,
              "Inconsistent tet_tri_neighbs size");

  ArgErrLogIf(tet_tet_neighbs.size() != pTetsN * 4,
              "Inconsistent tet_tet_neighbs size");

  // see comment in other constructor
  srand(time(nullptr));

  // copy vertex information and update bounding box
  pVerts.resize(pVertsN);
  for (size_t i = 0, j = 0; i < pVertsN; ++i, j += 3) {
    point3d vert = {verts[j], verts[j + 1], verts[j + 2]};

    pVerts[i] = vert;
    pBBox.insert(vert);
  }

  // copy tri information and compute tri barycenters.
  pTris.resize(pTrisN);
  pTri_norms.resize(pTrisN);
  pTri_barycs.resize(pTrisN);
  for (size_t i = 0, j = 0; i < pTrisN; ++i, j += 3) {
    tri_verts tri = {{tris[j], tris[j + 1], tris[j + 2]}};

    pTris[i] = tri;
    pTri_norms[i] = point3d{tri_norms[j], tri_norms[j + 1], tri_norms[j + 2]};
    pTri_barycs[i] = steps::math::tri_barycenter(
        pVerts[tri[0].get()], pVerts[tri[1].get()], pVerts[tri[2].get()]);
  }

  pTri_areas = tri_areas;
  for (auto a : pTri_areas) {

    ArgErrLogIf(a <= 0, "triangle with non-positive area");
  }
  pTri_diffboundaries.assign(pTrisN, nullptr);
  pTri_patches.assign(pTrisN, nullptr);

  // use tri data to make pBars, pTriBars.
  buildBarData();

  // copy tet information and compute tet barycenters.
  pTets.resize(pTetsN);
  pTet_barycenters.resize(pTetsN);
  for (size_t i = 0, j = 0; i < pTetsN; ++i, j += 4) {
    tet_verts tet = {{tets[j], tets[j + 1], tets[j + 2], tets[j + 3]}};

    pTets[i] = tet;
    pTet_barycenters[i] =
        steps::math::tet_barycenter(pVerts[tet[0].get()], pVerts[tet[1].get()],
                                    pVerts[tet[2].get()], pVerts[tet[3].get()]);
  }
  pTet_vols = tet_vols;
  for (auto v : pTet_vols) {
    ArgErrLogIf(v <= 0, "tetrahedron with non-positive volume");
  }

  // copy tet/tri and tet/tet neighbourhood infoirmation.
  pTri_tet_neighbours.resize(pTrisN);
  for (size_t i = 0, j = 0; i < pTrisN; ++i, j += 2) {
    pTri_tet_neighbours[i] =
        tri_tets{{tri_tet_neighbs[j], tri_tet_neighbs[j + 1]}};
  }

  pTet_tri_neighbours.resize(pTetsN);
  pTet_tet_neighbours.resize(pTetsN);
  for (size_t i = 0, j = 0; i < pTetsN; ++i, j += 4) {
    pTet_tri_neighbours[i] =
        tet_tris{{tet_tri_neighbs[j], tet_tri_neighbs[j + 1],
                  tet_tri_neighbs[j + 2], tet_tri_neighbs[j + 3]}};
    pTet_tet_neighbours[i] =
        tet_tets{{tet_tet_neighbs[j], tet_tet_neighbs[j + 1],
                  tet_tet_neighbs[j + 2], tet_tet_neighbs[j + 3]}};
  }

  pTet_comps.assign(pTetsN, nullptr);
}

////////////////////////////////////////////////////////////////////////////////

Tetmesh::~Tetmesh() {
  for (auto &membs : pMembs) {
    delete membs.second;
  }
  for (auto &diffb : pDiffBoundaries) {
    delete diffb.second;
  }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getVertex(vertex_id_t vidx) const {

  ArgErrLogIf(vidx >= pVertsN, "Vertex index is out of range.");

  auto const &vertex = pVerts[vidx.get()];
  return {vertex.begin(), vertex.end()};
}

////////////////////////////////////////////////////////////////////////////////

std::vector<tetrahedron_id_t>
Tetmesh::getVertexTetNeighbs(vertex_id_t vidx) const {
  ArgErrLogIf(vidx >= pVertsN, "Vertex index is out of range.");

  std::vector<tetrahedron_id_t> tets;
  for (auto tidx = 0u; tidx < pTetsN; ++tidx) {
    const tet_verts &v = pTets[tidx];
    if (v[0] == vidx || v[1] == vidx || v[2] == vidx || v[3] == vidx)
      tets.push_back(tidx);
  }
  return tets;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getBar(bar_id_t bidx) const {

  ArgErrLogIf(bidx >= pBarsN, "Bar index is out of range.");

  auto const &bar = pBars[bidx.get()];
  return strong_type_to_value_type(bar.begin(), bar.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTri(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  auto const &tri = pTris[tidx.get()];
  return strong_type_to_value_type(tri.begin(), tri.end());
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getTriArea(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  return pTri_areas[tidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getTriBarycenter(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  return as_vector(pTri_barycs[tidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getTriNorm(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  return as_vector(pTri_norms[tidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

TmPatch *Tetmesh::getTriPatch(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  return pTri_patches[tidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::setTriPatch(triangle_id_t tidx, TmPatch *patch) {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  pTri_patches[tidx.get()] = patch;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::setTriDiffBoundary(triangle_id_t tidx, DiffBoundary *diffb) {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  pTri_diffboundaries[tidx.get()] = diffb;
}

////////////////////////////////////////////////////////////////////////////////

DiffBoundary *Tetmesh::getTriDiffBoundary(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  return pTri_diffboundaries[tidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::setBarSDiffBoundary(bar_id_t bidx, SDiffBoundary *sdiffb) {

  ArgErrLogIf(bidx >= pBarsN, "Bar index is out of range.");

  pBar_sdiffboundaries[bidx.get()] = sdiffb;
}

////////////////////////////////////////////////////////////////////////////////

SDiffBoundary *Tetmesh::getBarSDiffBoundary(bar_id_t bidx) const {

  ArgErrLogIf(bidx >= pBarsN, "Bar index is out of range.");

  return pBar_sdiffboundaries[bidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTriBars(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  auto const &bars = pTri_bars[tidx.get()];
  return strong_type_to_value_type(bars.begin(), bars.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTriTetNeighb(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  auto const &neighbours = pTri_tet_neighbours[tidx.get()];
  return strong_type_to_value_type(neighbours.begin(), neighbours.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTriTriNeighb(triangle_id_t tidx,
                                              const TmPatch *tmpatch) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  // Triangles are neighbours if they share a bar
  std::vector<index_t> neighbours(3, triangle_id_t::unknown_value());
  tri_bars bars = pTri_bars[tidx.get()];

  for (triangle_id_t tri = 0u; tri < pTrisN; ++tri) {
    if (tri == tidx || pTri_patches[tri.get()] != tmpatch) {
      continue;
    }

    tri_bars neighbtribars = pTri_bars[tri.get()];
    int next = 0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (neighbtribars[j] != bars[i]) {
          continue;
        }

        if (neighbours[i] != triangle_id_t::unknown_value()) {
          std::ostringstream os;
          os << "Error in Patch initialisation for '" << tmpatch->getID()
             << "'. Patch triangle idx " << tidx
             << " found to have more than 3 neighbours.";
          ArgErrLog(os.str());
        }

        neighbours[i] = tri.get();
        next = 1;
        break;
      }
      if (next) {
        break;
      }
    }
  }
  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<index_t> Tetmesh::getTriTriNeighb(triangle_id_t tidx) const {

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  // Triangles are neighbours if they share a bar
  std::vector<index_t> neighbours(3, triangle_id_t::unknown_value());
  tri_bars bars = pTri_bars[tidx.get()];

  for (triangle_id_t tri = 0u; tri < pTrisN; ++tri) {
    if (tri == tidx || pTri_patches[tri.get()] == nullptr) {
      continue;
    }

    tri_bars neighbtribars = pTri_bars[tri.get()];
    int next = 0;
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (neighbtribars[j] != bars[i]) {
          continue;
        }

        if (neighbours[i] != triangle_id_t::unknown_value()) {
          std::ostringstream os;
          os << "Error: Triangle idx " << tidx
             << " found to have more than 3 neighbours.";
          ArgErrLog(os.str());
        }

        neighbours[i] = tri.get();
        next = 1;
        break;
      }
      if (next) {
        break;
      }
    }
  }
  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

std::set<index_t> Tetmesh::getTriTriNeighbs(triangle_id_t tidx) const {
  std::set<index_t> neighbours;

  ArgErrLogIf(tidx >= pTrisN, "Triangle index is out of range.");

  // Triangles are neighbours if they share a bar

  const tri_bars &bars = pTri_bars[tidx.get()];

  for (triangle_id_t tri = 0u; tri < pTrisN; ++tri) {
    if (tri == tidx) {
      continue;
    }

    const tri_bars &neighbtribars = pTri_bars[tri.get()];
    for (int i = 0; i < 3; ++i) {
      for (int j = 0; j < 3; ++j) {
        if (neighbtribars[j] != bars[i]) {
          continue;
        }

        neighbours.insert(tri.get());
        goto next_tri;
      }
    }
  next_tri:;
  }

  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

std::set<triangle_id_t> Tetmesh::getBarTriNeighbs(bar_id_t bidx) const {

  ArgErrLogIf(bidx >= pBarsN, "Bar index is out of range.");

  std::set<triangle_id_t> neighbours;

  for (uint tri = 0; tri < pTrisN; ++tri) {
    tri_bars neighbtribars = pTri_bars[tri];
    for (int i = 0; i < 3; ++i) {
      if (neighbtribars[i] != bidx) {
        continue;
      }

      neighbours.insert(tri);
      goto next_tri;
    }
  next_tri:;
  }

  return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::setBarTris(bar_id_t bidx, triangle_id_t itriidx,
                         triangle_id_t otriidx) {

  ArgErrLogIf(bidx >= pBarsN, "Bar index is out of range.");

  ArgErrLogIf(itriidx.unknown() or otriidx.unknown() or
                  itriidx >= pTrisN or otriidx >= pTrisN,
              "Invalid triangle index.");

  if (pBar_tri_neighbours[bidx.get()][0].valid() or pBar_tri_neighbours[bidx.get()][1].valid()) {
      std::ostringstream os;
      os << "Bar " << bidx << " is part of more than one surface diffusion boundary.";
      ArgErrLog(os.str());
  }

  pBar_tri_neighbours[bidx.get()][0] = itriidx;
  pBar_tri_neighbours[bidx.get()][1] = otriidx;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_flipTriTetNeighb(triangle_id_t tidx) {
  AssertLog(tidx < pTrisN);

  tri_tets &tt = pTri_tet_neighbours[tidx.get()];
  std::swap(tt[0], tt[1]);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_flipTriVerts(triangle_id_t tidx) {
  AssertLog(tidx < pTrisN);

  tri_verts &tri = pTris[tidx.get()];
  std::swap(tri[0], tri[1]);

  // also recalculate the triangle's normal; the long way (though should be
  // negative of previous normal)

  pTri_norms[tidx.get()] = steps::math::tri_normal(
      pVerts[tri[0].get()], pVerts[tri[1].get()], pVerts[tri[2].get()]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTet(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  auto const &tet = pTets[tidx.get()];
  return strong_type_to_value_type(tet.begin(), tet.end());
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getMeshVolume() const {
  return std::accumulate(pTet_vols.begin(), pTet_vols.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getTetVol(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  return pTet_vols[tidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getTetQualityRER(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  tet_verts verts = pTets[tidx.get()];
  const point3d &v0 = pVerts[verts[0].get()], &v1 = pVerts[verts[1].get()],
                &v2 = pVerts[verts[2].get()], &v3 = pVerts[verts[3].get()];

  return steps::math::tet_circumrad(v0, v1, v2, v3) /
         steps::math::tet_shortestedge(v0, v1, v2, v3);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getTetBarycenter(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  return as_vector(pTet_barycenters[tidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

TmComp *Tetmesh::getTetComp(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  return pTet_comps[tidx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::setTetComp(tetrahedron_id_t tidx, TmComp *comp) {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  pTet_comps[tidx.get()] = comp;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTetTriNeighb(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  auto const &tris = pTet_tri_neighbours[tidx.get()];
  return strong_type_to_value_type(tris.begin(), tris.end());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getTetTetNeighb(tetrahedron_id_t tidx) const {

  ArgErrLogIf(tidx >= pTetsN, "Tetrahedron index is out of range.");

  auto const &tets = pTet_tet_neighbours[tidx.get()];
  return strong_type_to_value_type(tets.begin(), tets.end());
}

////////////////////////////////////////////////////////////////////////////////

tetrahedron_id_t Tetmesh::findTetByPoint(std::vector<double> const &p) const {
  point3d x{p[0], p[1], p[2]};
  return findTetByPoint(x);
}

tetrahedron_id_t Tetmesh::findTetByPoint(point3d const &p) const {
  if (!pBBox.contains(p)) {
    return std::nullopt;
  }

  for (auto tidx = 0u; tidx < pTetsN; ++tidx) {
    const tet_verts &v = pTets[tidx];
    if (steps::math::tet_inside(pVerts[v[0].get()], pVerts[v[1].get()],
                                pVerts[v[2].get()], pVerts[v[3].get()], p)) {
      return tidx;
    }
  }

  return std::nullopt;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getBoundMin() const {
  return as_vector(pBBox.min());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getBoundMax() const {
  return as_vector(pBBox.max());
}

////////////////////////////////////////////////////////////////////////////////

// Created by weiliang 2010.02.02
std::vector<index_t> Tetmesh::getSurfTris() const {
  std::vector<index_t> tribounds;
  for (auto t = 0u; t < pTrisN; t++) {
    const tri_tets &trineighbor = pTri_tet_neighbours[t];
    if (trineighbor[0].unknown() || trineighbor[1].unknown()) {
      tribounds.push_back(t);
    }
  }
  return tribounds;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_checkMembID(std::string const &id) const {
  checkID(id);

  ArgErrLogIf(pMembs.find(id) != pMembs.end(),
              "'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleMembIDChange(std::string const &o, std::string const &n) {
  MembPMapCI m_old = pMembs.find(o);
  AssertLog(m_old != pMembs.end());

  if (o == n) {
    return;
  }
  _checkMembID(n);

  Memb *m = m_old->second;
  AssertLog(m != nullptr);
  pMembs.erase(m->getID()); // or s_old->first
  pMembs.insert(MembPMap::value_type(n, m));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleMembAdd(Memb *memb) {
  AssertLog(memb->getContainer() == this);
  _checkMembID(memb->getID());
  pMembs.insert(MembPMap::value_type(memb->getID(), memb));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleMembDel(Memb *memb) {
  AssertLog(memb->getContainer() == this);
  pMembs.erase(memb->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_checkDiffBoundaryID(std::string const &id) const {
  checkID(id);

  ArgErrLogIf(pDiffBoundaries.find(id) != pDiffBoundaries.end(),
              "'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleDiffBoundaryIDChange(std::string const &o,
                                          std::string const &n) {
  DiffBoundaryPMapCI db_old = pDiffBoundaries.find(o);
  AssertLog(db_old != pDiffBoundaries.end());

  if (o == n) {
    return;
  }
  _checkDiffBoundaryID(n);

  DiffBoundary *db = db_old->second;
  AssertLog(db != nullptr);
  pDiffBoundaries.erase(db->getID()); // or s_old->first
  pDiffBoundaries.insert(DiffBoundaryPMap::value_type(n, db));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleDiffBoundaryAdd(DiffBoundary *diffb) {
  AssertLog(diffb->getContainer() == this);
  _checkDiffBoundaryID(diffb->getID());
  pDiffBoundaries.insert(DiffBoundaryPMap::value_type(diffb->getID(), diffb));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleDiffBoundaryDel(DiffBoundary *diffb) {
  AssertLog(diffb->getContainer() == this);
  pDiffBoundaries.erase(diffb->getID());
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_checkSDiffBoundaryID(std::string const &id) const {
  checkID(id);

  ArgErrLogIf(pSDiffBoundaries.find(id) != pSDiffBoundaries.end(),
              "'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleSDiffBoundaryIDChange(std::string const &o,
                                           std::string const &n) {
  SDiffBoundaryPMapCI db_old = pSDiffBoundaries.find(o);
  AssertLog(db_old != pSDiffBoundaries.end());

  if (o == n) {
    return;
  }
  _checkSDiffBoundaryID(n);

  SDiffBoundary *db = db_old->second;
  AssertLog(db != nullptr);
  pSDiffBoundaries.erase(db->getID()); // or s_old->first
  pSDiffBoundaries.insert(SDiffBoundaryPMap::value_type(n, db));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleSDiffBoundaryAdd(SDiffBoundary *sdiffb) {
  AssertLog(sdiffb->getContainer() == this);
  _checkSDiffBoundaryID(sdiffb->getID());
  pSDiffBoundaries.insert(
      SDiffBoundaryPMap::value_type(sdiffb->getID(), sdiffb));
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::_handleSDiffBoundaryDel(SDiffBoundary *sdiffb) {
  AssertLog(sdiffb->getContainer() == this);
  pSDiffBoundaries.erase(sdiffb->getID());
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::SDiffBoundary *Tetmesh::_getSDiffBoundary(uint gidx) const {
  AssertLog(gidx < pSDiffBoundaries.size());
  auto db_it = pSDiffBoundaries.begin();
  std::advance(db_it, gidx);
  return db_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::Memb *Tetmesh::_getMemb(uint gidx) const {
  AssertLog(gidx < pMembs.size());
  auto mb_it = pMembs.begin();
  std::advance(mb_it, gidx);
  return mb_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::DiffBoundary *Tetmesh::_getDiffBoundary(uint gidx) const {
  AssertLog(gidx < pDiffBoundaries.size());
  auto db_it = pDiffBoundaries.begin();
  std::advance(db_it, gidx);
  return db_it->second;
}

////////////////////////////////////////////////////////////////////////
// Batch Data Access
////////////////////////////////////////////////////////////////////////

// Helper utility functions for batch data access:

/* Copy the components of n objects of type T to out_iter, where the objects
 * are specified by indices given by idx_iter over the collection items.
 */
template <typename T, typename I, typename J>
void batch_copy_components_n(
    const std::vector<T> &items, I idx_iter, size_t n, J out_iter,
    typename std::enable_if<std::is_pointer<I>::value>::type * = 0) {
  typename std::remove_const<typename std::remove_pointer<I>::type>::type index;
  try {
    for (size_t i = 0; i < n; ++i) {
      index = *idx_iter++;
      const auto &item = items.at(deref_strongid(index));
      out_iter = std::transform(
          item.begin(), item.end(), out_iter,
          [](const typename T::value_type &e) { return deref_strongid(e); });
    }
  } catch (const std::out_of_range &e) {
    ArgErrLog("Index out of range: no item with index " +
              std::to_string(index) + ".");
  }
}

template <typename T, typename I, typename J>
void batch_copy_components_n(
    const std::vector<T> &items, I idx_iter, size_t n, J out_iter,
    typename std::enable_if<!std::is_pointer<I>{}>::type * = 0) {
  typename std::iterator_traits<I>::value_type index;
  try {
    for (size_t i = 0; i < n; ++i) {
      index = *idx_iter++;
      const auto &item = items.at(deref_strongid(index));
      out_iter = std::transform(
          item.begin(), item.end(), out_iter,
          [](const typename T::value_type &e) { return deref_strongid(e); });
    }
  } catch (const std::out_of_range &e) {
    ArgErrLog("Index out of range: no item with index " +
              std::to_string(index) + ".");
  }
}

/* Copy n elements of the vector items to out_iter, where the elements
 * are specified by indices given by idx_iter.
 */

template <typename T, typename I, typename J>
void batch_copy_n(const std::vector<T> &items, I idx_iter, size_t n,
                  J out_iter) {
  size_t index = 0;
  try {
    for (size_t i = 0; i < n; ++i) {
      index = *idx_iter++;
      *out_iter++ = items.at(index);
    }
  } catch (const std::out_of_range &e) {
    ArgErrLog("Index out of range: no item with index " +
              std::to_string(index) + ".");
  }
}
////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getBatchTetBarycenters(
    std::vector<tetrahedron_id_t> const &tets) const {
  std::vector<double> data(tets.size() * 3);

  batch_copy_components_n(pTet_barycenters, tets.begin(), tets.size(),
                          data.begin());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTetBarycentersNP(const tetrahedron_id_t *indices,
                                       int input_size, double *centers,
                                       int output_size) const {

  ArgErrLogIf(input_size * 3 != output_size,
              "Length of output array should be 3 * length of input array.");

  batch_copy_components_n(pTet_barycenters, indices, input_size, centers);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Tetmesh::getBatchTriBarycenters(std::vector<triangle_id_t> const &tris) const {
  const auto ntris = tris.size();
  std::vector<double> data(ntris * 3);

  batch_copy_components_n(pTri_barycs, tris.begin(), ntris, data.begin());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTriBarycentersNP(const triangle_id_t *indices,
                                       int input_size, double *centers,
                                       int output_size) const {

  ArgErrLogIf(input_size * 3 != output_size,
              "Length of output array should be 3 * length of input array.");

  batch_copy_components_n(pTri_barycs, indices, input_size, centers);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Tetmesh::getBatchVertices(std::vector<index_t> const &verts) const {
  const auto nverts = verts.size();
  std::vector<double> data(nverts * 3);

  batch_copy_components_n(pVerts, verts.begin(), nverts, data.begin());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchVerticesNP(const index_t *indices, int input_size,
                                 double *coordinates, int output_size) const {

  ArgErrLogIf(input_size * 3 != output_size,
              "Length of output array (coordinates) should be 3 * length of "
              "input array (indices).");

  batch_copy_components_n(pVerts, indices, input_size, coordinates);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t>
Tetmesh::getBatchTris(std::vector<index_t> const &tris) const {
  std::vector<index_t> data(tris.size() * 3, 0);
  batch_copy_components_n(pTris, tris.begin(), tris.size(), data.begin());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTrisNP(const index_t *t_indices, int input_size,
                             index_t *v_indices, int output_size) const {

  ArgErrLogIf(input_size * 3 != output_size,
              "Length of output array should be 3 * length of input array.");

  batch_copy_components_n(pTris, t_indices, input_size, v_indices);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t>
Tetmesh::getBatchTets(std::vector<index_t> const &tets) const {
  const auto ntets = tets.size();
  std::vector<index_t> data(ntets * 4, 0);

  batch_copy_components_n(pTets, tets.begin(), ntets, data.begin());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTetsNP(const index_t *t_indices, int input_size,
                             index_t *v_indices, int output_size) const {

  ArgErrLogIf(input_size * 4 != output_size,
              "Length of output array should be 4 * length of input array.");

  batch_copy_components_n(pTets, t_indices, input_size, v_indices);
}

////////////////////////////////////////////////////////////////////////////////

uint Tetmesh::getTriVerticesSetSizeNP(const index_t *t_indices,
                                      int input_size) const {
  std::set<index_t> unique_indices;
  batch_copy_components_n(
      pTris, t_indices, input_size,
      std::inserter(unique_indices, unique_indices.begin()));

  return unique_indices.size();
}

////////////////////////////////////////////////////////////////////////////////

uint Tetmesh::getTetVerticesSetSizeNP(const index_t *t_indices,
                                      int input_size) const {
  std::set<tetrahedron_id_t ::value_type> unique_indices;
  batch_copy_components_n(
      pTets, t_indices, input_size,
      std::inserter(unique_indices, unique_indices.begin()));

  return unique_indices.size();
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getTriVerticesMappingSetNP(const index_t *t_indices,
                                         int input_size, index_t *t_vertices,
                                         int t_vertices_size, index_t *v_set,
                                         int v_set_size) const {

  ArgErrLogIf(
      input_size * 3 != t_vertices_size,
      "Length of t_vertices array should be 3 * length of input array.");

  // put unique vertices in v_set
  auto v_set_indices = make_unique_indexer<vertex_id_t>(v_set);
  int t_vertices_index = 0;

  for (int t = 0; t < input_size; ++t) {
    uint tidx = t_indices[t];

    ArgErrLogIf(tidx >= pTrisN, "Index out of range: no item with index " +
                                    std::to_string(tidx) + ".");

    for (auto vidx : pTris[tidx]) {
      t_vertices[t_vertices_index++] = v_set_indices[vidx];
    }
  }

  AssertLog(static_cast<int>(v_set_indices.size()) == v_set_size);
  AssertLog(t_vertices_index == t_vertices_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getTetVerticesMappingSetNP(const index_t *t_indices,
                                         int input_size, index_t *t_vertices,
                                         int t_vertices_size, index_t *v_set,
                                         int v_set_size) const {

  ArgErrLogIf(
      input_size * 4 != t_vertices_size,
      "Length of t_vertices array should be 4 * length of input array.");

  auto v_set_indices = make_unique_indexer<index_t>(v_set);
  int t_vertices_index = 0;

  for (int t = 0; t < input_size; ++t) {
    auto tidx = t_indices[t];

    ArgErrLogIf(tidx >= pTetsN, "Index out of range: no item with index " +
                                    std::to_string(tidx) + ".");

    for (auto vidx : pTets[tidx]) {
      t_vertices[t_vertices_index++] = v_set_indices[vidx.get()];
    }
  }

  AssertLog(static_cast<int>(v_set_indices.size()) == v_set_size);
  AssertLog(t_vertices_index == t_vertices_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genPointsInTet(tetrahedron_id_t tidx, unsigned npnts,
                             double *coords, unsigned int coord_size) const {

  ArgErrLogIf(npnts * 3 != coord_size,
              "Coordinate array size should be 3 * npnts.");

  ArgErrLogIf(tidx >= pTetsN, "Index out of range: no tetrahedron with index " +
                                  std::to_string(tidx) + ".");

  auto const &tet = pTets[tidx.get()];
  point3d v[4] = {pVerts[tet[0].get()], pVerts[tet[1].get()],
                  pVerts[tet[2].get()], pVerts[tet[3].get()]};

  for (uint p = 0, j = 0; p < npnts; p++, j += 3) {
    auto s = static_cast<double>(rand()) / RAND_MAX;
    auto t = static_cast<double>(rand()) / RAND_MAX;
    auto u = static_cast<double>(rand()) / RAND_MAX;
    point3d x = steps::math::tet_ranpnt(v[0], v[1], v[2], v[3], s, t, u);
    coords[j] = x[0];
    coords[j + 1] = x[1];
    coords[j + 2] = x[2];
  }
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genPointsInTri(triangle_id_t tidx, unsigned npnts, double *coords,
                             unsigned int coord_size) const {
  ArgErrLogIf(npnts * 3 != coord_size,
              "Coordinate array size should be 3 * npnts.");

  ArgErrLogIf(tidx >= pTrisN, "Index out of range: no triangle with index " +
                                  std::to_string(tidx) + ".");

  auto const &tri = pTris[tidx.get()];
  const point3d *v[3] = {&pVerts[tri[0].get()], &pVerts[tri[1].get()],
                         &pVerts[tri[2].get()]};

  for (uint p = 0, j = 0; p < npnts; p++, j += 3) {
    auto s = static_cast<double>(rand()) / RAND_MAX;
    auto t = static_cast<double>(rand()) / RAND_MAX;

    point3d x = steps::math::tri_ranpnt(*v[0], *v[1], *v[2], s, t);
    coords[j] = x[0];
    coords[j + 1] = x[1];
    coords[j + 2] = x[2];
  }
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genTetVisualPointsNP(const index_t *indices,
                                   unsigned int index_size,
                                   const unsigned int *point_counts,
                                   unsigned int count_size, double *coords,
                                   unsigned int coord_size) const {
  if (index_size != count_size) {
    ArgErrLogIf(
        index_size != count_size,
        "Length of point_counts array should be length of indices array.");
  }

  uint cpos = 0;
  for (auto i = 0u; i < index_size; i++) {
    auto npnts = point_counts[i];
    auto ncoords = 3 * npnts;
    ArgErrLogIf(cpos + ncoords > coord_size,
                "Length of coords array too short.");

    genPointsInTet(indices[i], npnts, coords + cpos, ncoords);
    cpos += ncoords;
  }

  ArgErrLogIf(cpos != coord_size,
              "Length of coords array longer than expected.");
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genTriVisualPointsNP(const index_t *indices,
                                   unsigned int index_size,
                                   const unsigned int *point_counts,
                                   unsigned int count_size, double *coords,
                                   unsigned int coord_size) const {

  ArgErrLogIf(
      index_size != count_size,
      "Length of point_counts array should be length of indices array.");

  uint cpos = 0;
  for (uint i = 0; i < index_size; i++) {
    uint npnts = point_counts[i];
    uint ncoords = 3 * npnts;

    ArgErrLogIf(cpos + ncoords > coord_size,
                "Length of coords array too short.");

    genPointsInTri(indices[i], npnts, coords + cpos, ncoords);
    cpos += ncoords;
  }

  ArgErrLogIf(cpos != coord_size,
              "Length of coords array longer than expected.");
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTetVolsNP(const index_t *indices, int index_size,
                                double *volumes, int volume_size) const {
  ArgErrLogIf(index_size != volume_size,
              "Length of volumes array should be length of indices array.");

  batch_copy_n(pTet_vols, indices, index_size, volumes);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getBatchTriAreasNP(const index_t *indices, int index_size,
                                 double *areas, int area_size) const {
  ArgErrLogIf(index_size != area_size,
              "Length of areas array should be length of indices array.");

  batch_copy_n(pTri_areas, indices, index_size, areas);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::reduceBatchTetPointCountsNP(const index_t *indices,
                                          unsigned int index_size,
                                          unsigned int *point_counts,
                                          unsigned int count_size,
                                          double max_density) {

  ArgErrLogIf(
      index_size != count_size,
      "Length of point_counts array should be length of indices array.");

  for (auto i = 0u; i < index_size; i++) {
    auto tidx = indices[i];

    ArgErrLogIf(tidx >= pTetsN,
                "Index out of range: no tetrahedron with index " +
                    std::to_string(tidx) + ".");

    point_counts[i] = std::min(
        point_counts[i], static_cast<uint>(max_density * pTet_vols[tidx]));
  }
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::reduceBatchTriPointCountsNP(const index_t *indices,
                                          unsigned int index_size,
                                          unsigned int *point_counts,
                                          unsigned int count_size,
                                          double max_density) {

  ArgErrLogIf(
      index_size != count_size,
      "Length of point_counts array should be length of indices array.");

  for (auto i = 0u; i < index_size; i++) {
    uint tidx = indices[i];
    ArgErrLogIf(tidx >= pTrisN, "Index out of range: no triangle with index " +
                                    std::to_string(tidx) + ".");

    point_counts[i] = std::min(
        point_counts[i], static_cast<uint>(max_density * pTri_areas[tidx]));
  }
}

////////////////////////////////////////////////////////////////////////
// ROI (Region of Interest) Data
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::addROI(std::string const &id, steps::tetmesh::ElementType type,
                     const ROISet::set_data_type &indices) {
  bool inserted = false;
  switch (type) {
  case ELEM_VERTEX: {
    ROITypeTraits<ROI_VERTEX>::data_type data(indices.begin(), indices.end());
    inserted = rois.insert<ROI_VERTEX>(id, data).second;
    break;
  }
  case ELEM_TRI: {
    ROITypeTraits<ROI_TRI>::data_type data(indices.begin(), indices.end());
    inserted = rois.insert<ROI_TRI>(id, data).second;
    break;
  }
  case ELEM_TET: {
    ROITypeTraits<ROI_TET>::data_type data(indices.begin(), indices.end());
    inserted = rois.insert<ROI_TET>(id, data).second;
    break;
  }
  case ELEM_UNDEFINED:
    break;
  }
  if (!inserted) {
    CLOG(WARNING, "general_log")
        << "ROI data with id " << id
        << " already exists. Use replaceROI() to replace the data.\n";
  }
}

////////////////////////////////////////////////////////////////////////////////
void Tetmesh::removeROI(std::string const &id) {
  if (rois.erase<ROI_TRI>(id) != 0) {
    return;
  }
  if (rois.erase<ROI_TET>(id) != 0) {
    return;
  }
  if (rois.erase<ROI_VERTEX>(id) != 0) {
    return;
  }
  CLOG(WARNING, "general_log")
      << "Unable to find ROI data with id " << id << ".\n";
}

////////////////////////////////////////////////////////////////////////////////
void Tetmesh::replaceROI(std::string const &id,
                         steps::tetmesh::ElementType type,
                         const ROISet::set_data_type &indices) {
  bool replaced = false;
  switch (type) {
  case ELEM_VERTEX: {
    ROITypeTraits<ROI_VERTEX>::data_type data(indices.begin(), indices.end());
    replaced = rois.replace<ROI_VERTEX>(id, data);
    break;
  }
  case ELEM_TRI: {
    ROITypeTraits<ROI_TRI>::data_type data(indices.begin(), indices.end());
    replaced = rois.replace<ROI_TRI>(id, data);
    break;
  }
  case ELEM_TET: {
    ROITypeTraits<ROI_TET>::data_type data(indices.begin(), indices.end());
    replaced = rois.replace<ROI_TET>(id, data);
    break;
  }
  case ELEM_UNDEFINED:
    break;
  }
  if (!replaced) {
    CLOG(WARNING, "general_log")
        << "Unable to find ROI data with id " << id << ".\n";
  }
}

////////////////////////////////////////////////////////////////////////////////
steps::tetmesh::ElementType Tetmesh::getROIType(std::string const &id) const {
  if (rois.get<ROI_VERTEX>(id, 0, false) != rois.end<ROI_VERTEX>()) {
    return ELEM_VERTEX;
  }
  if (rois.get<ROI_TRI>(id, 0, false) != rois.end<ROI_TRI>()) {
    return ELEM_TRI;
  }
  if (rois.get<ROI_TET>(id, 0, false) != rois.end<ROI_TET>()) {
    return ELEM_TET;
  }

  CLOG(WARNING, "general_log")
      << "Unable to find ROI data with id " << id << ".\n";
  return steps::tetmesh::ELEM_UNDEFINED;
}

////////////////////////////////////////////////////////////////////////////////
ROISet::vector_data_type Tetmesh::getROIData(std::string const &id) const {
  {
    auto const &roi = rois.get<ROI_TRI>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TRI>()) {
      return strong_type_to_value_type(roi->second.begin(), roi->second.end());
    }
  }
  {
    auto const &roi = rois.get<ROI_TET>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TET>()) {
      return strong_type_to_value_type(roi->second.begin(), roi->second.end());
    }
  }
  {
    auto const &roi =
        rois.get<ROI_VERTEX>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_VERTEX>()) {
      return strong_type_to_value_type(roi->second.begin(), roi->second.end());
    }
  }
  CLOG(WARNING, "general_log")
      << "Unable to find ROI data with id " << id << ".\n";
  static const ROISet::vector_data_type empty;
  return empty;
}

////////////////////////////////////////////////////////////////////////////////
uint Tetmesh::getROIDataSize(std::string const &id) const {
  {
    auto const &roi = rois.get<ROI_TRI>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TRI>()) {
      return roi->second.size();
    }
  }
  {
    auto const &roi = rois.get<ROI_TET>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TET>()) {
      return roi->second.size();
    }
  }
  {
    auto const &roi =
        rois.get<ROI_VERTEX>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_VERTEX>()) {
      return roi->second.size();
    }
  }

  CLOG(WARNING, "general_log")
      << "Unable to find ROI data with id " << id << ".\n";
  return 0;
}

////////////////////////////////////////////////////////////////////////////////

uint Tetmesh::getNROIs() const { return rois.size(); }

////////////////////////////////////////////////////////////////////////////////

ROISet Tetmesh::getROI(std::string const &id) const {
  ROISet eax;
  {
    auto const &roi = rois.get<ROI_TRI>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TRI>()) {
      eax.type = ELEM_TRI;
      eax.indices.resize(roi->second.size());
      std::transform(roi->second.begin(), roi->second.end(),
                     eax.indices.begin(), deref_strong_id_func{});
      return eax;
    }
  }
  {
    auto const &roi = rois.get<ROI_TET>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_TET>()) {
      eax.type = ELEM_TET;
      eax.indices.resize(roi->second.size());
      std::transform(roi->second.begin(), roi->second.end(),
                     eax.indices.begin(), deref_strong_id_func{});
      return eax;
    }
  }
  {
    auto const &roi =
        rois.get<ROI_VERTEX>(id, 0 /* count */, false /* warning */);
    if (roi != rois.end<ROI_VERTEX>()) {
      eax.type = ELEM_VERTEX;
      eax.indices.resize(roi->second.size());
      std::transform(roi->second.begin(), roi->second.end(),
                     eax.indices.begin(), deref_strong_id_func{});
      return eax;
    }
  }
  CLOG(WARNING, "general_log")
      << "Unable to find ROI data with id " << id << ".\n";
  static const ROISet empty;
  return empty;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> Tetmesh::getAllROINames() const {
  std::vector<std::string> output;
  rois.ids(output);
  return output;
}

////////////////////////////////////////////////////////////////////////////////

bool Tetmesh::checkROI(std::string const &id, ElementType type, uint count,
                       bool warning) const {
  switch (type) {
  case ELEM_VERTEX:
    return rois.get<ROI_VERTEX>(id, count, warning) != rois.end<ROI_VERTEX>();
  case ELEM_TET:
    return rois.get<ROI_TET>(id, count, warning) != rois.end<ROI_TET>();
  case ELEM_TRI:
    return rois.get<ROI_TRI>(id, count, warning) != rois.end<ROI_TRI>();
  default:
    return false;
  }
}

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Tetmesh::getROITetBarycenters(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  const auto size = roi->second.size();
  std::vector<double> data(size * 3, 0.0);
  getBatchTetBarycentersNP(roi->second.data(), size, &data.front(),
                           data.size());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITetBarycentersNP(std::string const &ROI_id, double *centers,
                                     int output_size) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id, output_size / 3);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTetBarycentersNP(roi->second.data(), roi->second.size(), centers,
                           output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double>
Tetmesh::getROITriBarycenters(std::string const &ROI_id) const {
  const auto &roi = rois.get<ROI_TRI>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  const auto size = roi->second.size();
  std::vector<double> data(size * 3, 0.0);
  getBatchTriBarycentersNP(roi->second.data(), size, &data.front(),
                           data.size());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITriBarycentersNP(std::string const &ROI_id, double *centers,
                                     int output_size) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id, output_size / 3);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTriBarycentersNP(roi->second.data(), roi->second.size(), centers,
                           output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> Tetmesh::getROIVertices(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_VERTEX>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_VERTEX>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  const auto size = roi->second.size();
  std::vector<double> data(size * 3, 0.0);
  getBatchVerticesNP(reinterpret_cast<const index_t *>(roi->second.data()),
                     size, &data.front(), data.size());

  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROIVerticesNP(std::string const &ROI_id, double *coordinates,
                               int output_size) const {
  auto const &roi = rois.get<ROI_VERTEX>(ROI_id, output_size / 3);

  ArgErrLogIf(
      roi == rois.end<ROI_VERTEX>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchVerticesNP(reinterpret_cast<const index_t *>(roi->second.data()),
                     roi->second.size(), coordinates, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getROITris(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  const auto size = roi->second.size();
  std::vector<index_t> data(size * 3, 0.0);
  getBatchTrisNP(reinterpret_cast<const index_t *>(roi->second.data()), size,
                 &data.front(), data.size());

  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITrisNP(std::string const &ROI_id, index_t *v_indices,
                           int output_size) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTrisNP(reinterpret_cast<const index_t *>(roi->second.data()),
                 roi->second.size(), v_indices, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<index_t> Tetmesh::getROITets(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  const auto size = roi->second.size();
  std::vector<index_t> data(size * 4, 0.0);
  getBatchTetsNP(reinterpret_cast<const index_t *>(roi->second.data()), size,
                 &data.front(), data.size());
  return data;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITetsNP(std::string const &ROI_id, index_t *v_indices,
                           int output_size) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTetsNP(reinterpret_cast<const index_t *>(roi->second.data()),
                 roi->second.size(), v_indices, output_size);
}

////////////////////////////////////////////////////////////////////////////////

uint Tetmesh::getROITriVerticesSetSizeNP(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  return getTriVerticesSetSizeNP(
      reinterpret_cast<const index_t *>(roi->second.data()),
      roi->second.size());
}

////////////////////////////////////////////////////////////////////////////////

uint Tetmesh::getROITetVerticesSetSizeNP(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  return getTetVerticesSetSizeNP(
      reinterpret_cast<const index_t *>(roi->second.data()),
      roi->second.size());
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITriVerticesMappingSetNP(std::string const &ROI_id,
                                            index_t *t_vertices,
                                            int t_vertices_size, index_t *v_set,
                                            int v_set_size) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id, t_vertices_size / 3);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getTriVerticesMappingSetNP(
      reinterpret_cast<const index_t *>(roi->second.data()), roi->second.size(),
      t_vertices, t_vertices_size, v_set, v_set_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITetVerticesMappingSetNP(std::string const &ROI_id,
                                            index_t *t_vertices,
                                            int t_vertices_size, index_t *v_set,
                                            int v_set_size) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id, t_vertices_size / 4);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getTetVerticesMappingSetNP(
      reinterpret_cast<const index_t *>(roi->second.data()), roi->second.size(),
      t_vertices, t_vertices_size, v_set, v_set_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genROITetVisualPointsNP(std::string const &ROI_id,
                                      unsigned int *point_counts,
                                      int count_size, double *coords,
                                      int coord_size) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id, count_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  genTetVisualPointsNP(reinterpret_cast<const index_t *>(roi->second.data()),
                       roi->second.size(), point_counts, count_size, coords,
                       coord_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::genROITriVisualPointsNP(std::string const &ROI_id,
                                      unsigned int *point_counts,
                                      int count_size, double *coords,
                                      int coord_size) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id, count_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  genTriVisualPointsNP(reinterpret_cast<const index_t *>(roi->second.data()),
                       roi->second.size(), point_counts, count_size, coords,
                       coord_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITetVolsNP(std::string const &ROI_id, double *volumes,
                              int volume_size) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id, volume_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTetVolsNP(reinterpret_cast<const index_t *>(roi->second.data()),
                    roi->second.size(), volumes, volume_size);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::getROITriAreasNP(std::string const &ROI_id, double *areas,
                               int area_size) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id, area_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  getBatchTriAreasNP(reinterpret_cast<const index_t *>(roi->second.data()),
                     roi->second.size(), areas, area_size);
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getROIVol(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TET>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  double sum_vol = 0.0;
  for (auto &tet : roi->second) {
    sum_vol += getTetVol(tet);
  }

  return sum_vol;
}

////////////////////////////////////////////////////////////////////////////////

double Tetmesh::getROIArea(std::string const &ROI_id) const {
  auto const &roi = rois.get<ROI_TRI>(ROI_id);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  double sum_area = 0.0;
  for (const auto &tidx : roi->second) {
    sum_area += getTriArea(tidx);
  }
  return sum_area;
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::reduceROITetPointCountsNP(std::string const &ROI_id,
                                        unsigned int *point_counts,
                                        int count_size, double max_density) {
  auto const &roi = rois.get<ROI_TET>(ROI_id, count_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TET>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  reduceBatchTetPointCountsNP(
      reinterpret_cast<const index_t *>(roi->second.data()), roi->second.size(),
      point_counts, count_size, max_density);
}

////////////////////////////////////////////////////////////////////////////////

void Tetmesh::reduceROITriPointCountsNP(std::string const &ROI_id,
                                        unsigned int *point_counts,
                                        int count_size, double max_density) {
  auto const &roi = rois.get<ROI_TRI>(ROI_id, count_size);

  ArgErrLogIf(
      roi == rois.end<ROI_TRI>(),
      "ROI check fail, please make sure the ROI stores correct elements.");

  reduceBatchTriPointCountsNP(
      reinterpret_cast<const index_t *>(roi->second.data()), roi->second.size(),
      point_counts, count_size, max_density);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<tetrahedron_id_t, double>>
Tetmesh::intersectMontecarlo(const point3d &p_start, const point3d &p_end,
                             const tetrahedron_id_t &tet_start,
                             unsigned int sampling) const {
  std::vector<std::pair<tetrahedron_id_t, double>> segment_intersects;

  tetrahedron_id_t cur_tet =
      (tet_start.unknown()) ? findTetByPoint(p_start) : tet_start;
  if (cur_tet.unknown()) {
    CLOG(WARNING, "general_log") << "Initial point is not in the mesh.\n";
    return segment_intersects;
  }

  const vertex_id_t *tet_vertices = _getTet(cur_tet);

  // Test the case where 2 points are inside the same tet
  {
    if (tet_inside(_getVertex(tet_vertices[0]), _getVertex(tet_vertices[1]),
                   _getVertex(tet_vertices[2]), _getVertex(tet_vertices[3]),
                   p_start) &&
        (tet_inside(_getVertex(tet_vertices[0]), _getVertex(tet_vertices[1]),
                    _getVertex(tet_vertices[2]), _getVertex(tet_vertices[3]),
                    p_end))) {
      return {std::make_pair<>(cur_tet, 1.0)};
    }
  }

  const point3d p_diff = p_end - p_start;

  std::unordered_map<tetrahedron_id_t, int> accu;

  for (unsigned n = 0; n < sampling; ++n) {
    const point3d p_cur = p_start + p_diff * n / (sampling - 1);
    // not sure you need to check for UNKNOWN_TET here
    // check if sample point is in current tet
    if (cur_tet.valid() and !tet_inside(_getVertex(tet_vertices[0]),
                                        _getVertex(tet_vertices[1]),
                                        _getVertex(tet_vertices[2]),
                                        _getVertex(tet_vertices[3]),
                                        p_cur)) {
        // if not loop over neighbouring tets
        auto neighbTets = _getTetTetNeighb(cur_tet);
        cur_tet = std::nullopt;
        for (unsigned i = 0; i < 4; ++i) {
            if (neighbTets[i].unknown()) {
                continue;
            }
            tet_vertices = _getTet(neighbTets[i]);
            if (tet_inside(_getVertex(tet_vertices[0]),
                           _getVertex(tet_vertices[1]),
                           _getVertex(tet_vertices[2]),
                           _getVertex(tet_vertices[3]),
                           p_cur)) {
                cur_tet = neighbTets[i];
                break;
            }
        }
    }

    // failed to find point within neighbours, search all
    if (cur_tet.unknown()) {
      cur_tet = findTetByPoint(p_cur);
    }

    // FAIL
    if (cur_tet.unknown()) {
      CLOG(WARNING, "general_log") << "Could not find sample point.\n";
      continue;
    }

    tet_vertices = _getTet(cur_tet);

    accu[cur_tet] = accu[cur_tet] + 1;
  }

  for (const auto &kv : accu) {
    segment_intersects.emplace_back(kv.first,
                                    static_cast<double>(kv.second) / sampling);
  }

  return segment_intersects;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<std::pair<tetrahedron_id_t, double>>
Tetmesh::intersectDeterministic(const point3d &p_start, const point3d &p_end,
                                const tetrahedron_id_t &tet_start) const {
  std::vector<std::pair<tetrahedron_id_t, double>> segment_intersects;

  tetrahedron_id_t cur_tet =
      (tet_start.unknown()) ? findTetByPoint(p_start) : tet_start;
  if (cur_tet.unknown()) {
    CLOG(WARNING, "general_log") << "Initial point is not in the mesh.\n";
    return segment_intersects;
  }

  double seg_len = p_start.distance(p_end);
  point3d previous_point = p_start;
  const vertex_id_t *tet_vertices = _getTet(cur_tet);

  // storage for tet_neighbors
  std::vector<tetrahedron_id_t> tet_neighbors;

  // Loop over neighbor tets until we reach the tet containing segment end (or
  // outside) Will not even enter if the cur_tet contains p_end
  size_t neighbs_id = 0;
  while (!math::tet_inside(
      _getVertex(tet_vertices[0]), _getVertex(tet_vertices[1]),
      _getVertex(tet_vertices[2]), _getVertex(tet_vertices[3]), p_end)) {
    point3d intersection;
    tetrahedron_id_t next_tet = cur_tet;

    // Loop over current tet faces
    // if intersection is found select next tet
    for (const auto &tri_id : pTet_tri_neighbours[cur_tet.get()]) {
      const auto tri_v_ids = _getTri(tri_id);
      if (math::tri_intersect_line(
              _getVertex(tri_v_ids[0]), _getVertex(tri_v_ids[1]),
              _getVertex(tri_v_ids[2]), previous_point, p_end, intersection)) {

        // Check that intersection is new
        if (intersection.almostEqual(previous_point))
          continue;

        // Next tet is the other side of the chosen triangle
        auto tet_ids = _getTriTetNeighb(tri_id);
        next_tet = (tet_ids[0] == cur_tet) ? tet_ids[1] : tet_ids[0];
        break;
      }
    }

    // If next_tet == cur_tet then the line did not cross any face
    // this can happen due to finite precision arithmetic
    // try other neighboring tets
    if (next_tet == cur_tet) {
      // First pass, need to build neighbors list
      if (neighbs_id == 0) {
        for (size_t v = 0; v < 4; ++v) {
          // Get all tets sharing this vertex
          auto tet_v = getVertexTetNeighbs(tet_vertices[v]);
          // Store all neighbours except cur_tet
          for (size_t t = 0; t < tet_v.size(); ++t) {
            if (tet_v[t] != cur_tet)
              tet_neighbors.push_back(tet_v[t]);
          }
        }
      }

      // We tried all neighbors already, give up
      if (neighbs_id >= tet_neighbors.size()) {
        CLOG(WARNING, "general_log") << "Could not find endpoint.\n";
        return segment_intersects;
      }

      next_tet = tet_neighbors[neighbs_id++];
    }

    double ratio = previous_point.distance(intersection) / seg_len;
    if (ratio > math::tol_lin)
      segment_intersects.emplace_back(cur_tet, ratio);

    // Otherwise, unknown() then the segment finishes outside the mesh
    if (next_tet.unknown()) {
      CLOG(WARNING, "general_log") << "Endpoint is outside mesh.\n";
      return segment_intersects;
    }

    previous_point = intersection;
    cur_tet = next_tet;
    tet_vertices = _getTet(cur_tet);
  }

  // last intersection to end
  if (p_start == previous_point) {
    segment_intersects.emplace_back(cur_tet, 1.0);
  } else {
    double ratio = previous_point.distance(p_end) / seg_len;
    segment_intersects.emplace_back(cur_tet, ratio);
  }

  return segment_intersects;
}

std::vector<Tetmesh::intersection_list_t>
Tetmesh::intersect(const double *points, int n_points, int sampling) const {
  std::vector<Tetmesh::intersection_list_t> intersecs;
  if (n_points <= 1) {
    CLOG(WARNING, "general_log")
        << "Please provide at least two points to define a segment.\n";
    return intersecs;
  }

  point3d start_p(points[0], points[1],
                  points[2]); // beginning of the first segment
  tetrahedron_id_t cur_tet = findTetByPoint(start_p);

  // loop over each segment
  for (int i = 1; i < n_points; i++) {
    point3d end_p(points[i * 3], points[i * 3 + 1], points[i * 3 + 2]);
    std::vector<std::pair<tetrahedron_id_t, double>> isecs;
    if (sampling > 0)
      isecs = intersectMontecarlo(start_p, end_p, cur_tet, sampling);
    else
      isecs = intersectDeterministic(start_p, end_p, cur_tet);
    // Get last_tet right away. isecs gets moved later
    cur_tet = isecs.empty() ? std::nullopt : isecs.back().first;
    // We reinterpret cast since the types are equivallent but not auto
    // convertible
    intersecs.push_back(
        std::move(*reinterpret_cast<Tetmesh::intersection_list_t *>(&isecs)));
    start_p = end_p;
  }
  return intersecs;
}

} // namespace tetmesh
} // namespace steps
