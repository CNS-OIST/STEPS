/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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

#include <algorithm>
#include <cassert>
#include <cmath>  
#include <cstdlib>
#include <ctime>
#include <iostream>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <vector>
#include <numeric>

#include "steps/common.h"
#include "steps/error.hpp"

#include "steps/math/bbox.hpp"
#include "steps/math/point.hpp"
#include "steps/math/smallsort.hpp"
#include "steps/math/tetrahedron.hpp"
#include "steps/math/triangle.hpp"

#include "steps/geom/tetmesh.hpp"
#include "steps/geom/memb.hpp"

#include "steps/util/collections.hpp"
#include "steps/util/checkid.hpp"

namespace stetmesh = steps::tetmesh;

using steps::util::as_vector;
using steps::util::make_unique_indexer;
using steps::util::checkID;

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::Tetmesh(std::vector<double> const & verts,
                           std::vector<uint> const & tets,
                           std::vector<uint> const & tris)
: Geom()
, pVertsN(0)
, pBarsN(0)
, pTrisN(0)
, pTetsN(0)
, pMembs()
, pDiffBoundaries()
, pBar_tri_neighbours()
{
    using steps::math::small_sort;
    using std::to_string;

    srand (time(NULL));
    
    // check the vectors are of the expected size
    if (verts.size() % 3 || tris.size() % 3 || tets.size() % 4)
        throw ArgErr("Tables supplied to Tet mesh initialiser function are not of the expected dimensions");

    pVertsN = verts.size() / 3;
    if (!pVertsN) throw ArgErr("Empty vertex list");

    pTetsN = tets.size() / 4;
    if (!pTetsN) throw ArgErr("Empty tets list");

    // copy vertices and update bounding box
    pVerts.resize(pVertsN);
    for (uint i = 0, j = 0; i < pVertsN; ++i, j+=3) {
        point3d vert = {verts[j],verts[j+1],verts[j+2]};
        pVerts[i] = vert;
        pBBox.insert(vert);
    }

    // copy tets; 
    pTets.resize(pTetsN);
    for (uint i = 0, j = 0; i < pTetsN; ++i, j+=4)
        pTets[i] = tet_verts{tets[j],tets[j+1],tets[j+2],tets[j+3]};

    // Add user-supplied tris and faces for each tet to pTris and set
    // tet->tri adjacency
    {
        pTet_tri_neighbours.resize(pTetsN);
        auto tri_indices = make_unique_indexer<tri_verts>(std::back_inserter(pTris));

        // first add user-supplied tris:
        size_t userTrisN = tris.size()/3;
        for (uint i = 0; i < tris.size(); i+=3)
            tri_indices.insert(small_sort<3>(tri_verts{tris[i], tris[i+1], tris[i+2]}));

        for (uint i = 0; i < pTetsN; ++i) {
            const tet_verts &tet=pTets[i];
            tri_verts tris[4] = {
                small_sort<3>(tri_verts{tet[0],tet[1],tet[2]}),
                small_sort<3>(tri_verts{tet[0],tet[1],tet[3]}),
                small_sort<3>(tri_verts{tet[0],tet[2],tet[3]}),
                small_sort<3>(tri_verts{tet[1],tet[2],tet[3]})
            };

            for (int j = 0; j < 4; ++j)
                pTet_tri_neighbours[i][j] = tri_indices[tris[j]];
        }
    }
    pTris.shrink_to_fit();
    pTrisN = pTris.size();

    // For each tet, compute volume and barycentre, and update
    // tri->tet adjacency and tet->tet adjacency information.

    pTet_vols.resize(pTetsN);
    pTet_barycentres.resize(pTetsN);

    pTet_tet_neighbours.assign(pTetsN,tet_tets{-1,-1,-1,-1});
    pTri_tet_neighbours.assign(pTrisN,tri_tets{-1,-1});

    for (int i = 0; i < pTetsN; ++i) {
        auto tet = pTets[i];
        point3d v[4] = {pVerts[tet[0]], pVerts[tet[1]], pVerts[tet[2]], pVerts[tet[3]]};

        pTet_vols[i] = steps::math::tet_vol(v[0],v[1],v[2],v[3]);
        pTet_barycentres[i] = steps::math::tet_barycenter(v[0],v[1],v[2],v[3]);

        if (pTet_vols[i]<=0) throw ArgErr("degenerate tetrahedron "+to_string(i));

        for (int face = 0; face < 4; ++face) {
            int tri = pTet_tri_neighbours[i][face];
            auto &tri_tets = pTri_tet_neighbours[tri];

            if (tri_tets[0] == -1) tri_tets[0]=i;
            else if (tri_tets[1] == -1) {
                tri_tets[1] = i;

                // tri_tets[0] and [1] are neighbours.
                int other_tet = tri_tets[0];
                pTet_tet_neighbours[i][face]=other_tet;

                bool found_other_face = false;
                for (int other_face = 0; other_face < 4; ++other_face) {
                    if (pTet_tri_neighbours[other_tet][other_face] == tri) {
                        if (pTet_tet_neighbours[other_tet][other_face] != -1)
                            throw std::logic_error("inconsistent tet<->tet association");

                        pTet_tet_neighbours[other_tet][other_face] = i;
                        found_other_face = true;
                        break;
                    }
                }
                if (!found_other_face) 
                    throw std::logic_error("inconsistent tet<->tet association");

            }
            else throw std::logic_error("inconsisent tri<->tet association");
        }
    }

    // for each tri, compute area, barycentre, norm; allocate vectors for patches and diff. boundaries.

    pTri_areas.resize(pTrisN);
    pTri_barycs.resize(pTrisN);
    pTri_norms.resize(pTrisN);

    pTri_diffboundaries.assign(pTrisN,nullptr);
    pTri_patches.assign(pTrisN,nullptr);

    for (int i = 0; i < pTrisN; ++i) {
        auto tri = pTris[i];
        point3d v[3] = {pVerts[tri[0]], pVerts[tri[1]], pVerts[tri[2]]};

        pTri_areas[i] = steps::math::tri_area(v[0], v[1], v[2]);
        pTri_norms[i] = steps::math::tri_normal(v[0], v[1], v[2]);
        pTri_barycs[i] = steps::math::tri_barycenter(v[0], v[1], v[2]);

        if (pTri_areas[i]<=0) throw ArgErr("degenerate triangle "+to_string(i));
    }

    // initialise tet compartment association.
    pTet_comps.assign(pTetsN,nullptr);

    // create tri->bar adjacency
    buildBarData();

}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::buildBarData() {
    using steps::math::small_sort;

    pTri_bars.resize(pTrisN);

    pBars.clear();
    auto bar_indices = make_unique_indexer<bar_verts>(std::back_inserter(pBars));

    for (uint i = 0; i < pTrisN; ++i) {
        const tri_verts &tri=pTris[i];
        bar_verts bars[3] = {
            small_sort<2>(bar_verts{tri[0],tri[1]}),
            small_sort<2>(bar_verts{tri[0],tri[2]}),
            small_sort<2>(bar_verts{tri[1],tri[2]}),
        };

        for (int j = 0; j < 3; ++j) 
            pTri_bars[i][j] = bar_indices[bars[j]];
    }
    pBars.shrink_to_fit();
    pBarsN = pBars.size();

    // Include surface diffusion boundary stuff here
    pBar_sdiffboundaries.assign(pBarsN,nullptr);

    pBar_tri_neighbours.resize(pBarsN);
    pBar_tri_neighbours.assign(pBarsN,bar_tris{-1,-1});

}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::Tetmesh(std::vector<double> const & verts,
        std::vector<uint> const & tris,
        std::vector<double> const & tri_areas,
        std::vector<double> const & tri_norms,
        std::vector<int> const & tri_tet_neighbs,
        std::vector<uint> const & tets,
        std::vector<double> const & tet_vols,
        std::vector<double> const & tet_barycs,
        std::vector<uint> const & tet_tri_neighbs,
        std::vector<int> const & tet_tet_neighbs)
: Geom()
, pVertsN()
, pBarsN(0)
, pTrisN(0)
, pTetsN(0)
, pMembs()
, pDiffBoundaries()
{
    // check the vectors are of the expected size
    if (verts.size() % 3 || tris.size() % 3 || tets.size() % 4)
        throw ArgErr("Tables supplied to Tet mesh initialiser function are not of the expected dimensions");

    pVertsN = verts.size() / 3;
    if (!pVertsN) throw ArgErr("Empty vertex list");

    pTetsN = tets.size() / 4;
    if (!pTetsN) throw ArgErr("Empty tets list");

    pTrisN = tris.size() /3;
    if (!pTetsN) throw ArgErr("Empty tris list");

    if (tri_areas.size() != pTrisN)         throw ArgErr("Inconsistent tri_areas size");
    if (tri_norms.size() != pTrisN*3)       throw ArgErr("Inconsistent tri_norms size");
    if (tri_tet_neighbs.size() != pTrisN*2) throw ArgErr("Inconsistent tri_tet_neighbrs size");
    if (tet_vols.size() != pTetsN)          throw ArgErr("Inconsistent tet_vols size");
    if (tet_barycs.size() != pTetsN*3)      throw ArgErr("Inconsistent tet_barycs size");
    if (tet_tri_neighbs.size() != pTetsN*4) throw ArgErr("Inconsistent tet_tri_neighbs size");
    if (tet_tet_neighbs.size() != pTetsN*4) throw ArgErr("Inconsistent tet_tet_neighbs size");


    // see comment in other constructor
    srand (time(NULL));
    
    // copy vertex information and update bounding box
    pVerts.resize(pVertsN);
    for (uint i = 0, j = 0; i < pVertsN; ++i, j += 3) {
        point3d vert = {verts[j], verts[j+1], verts[j+2]};

        pVerts[i] = vert;
        pBBox.insert(vert);
    }

    // copy tri information and compute tri barycentres.
    pTris.resize(pTrisN);
    pTri_norms.resize(pTrisN);
    pTri_barycs.resize(pTrisN);
    for (uint i = 0, j = 0; i < pTrisN; ++i, j += 3) {
        tri_verts tri = {tris[j], tris[j+1], tris[j+2]};

        pTris[i] = tri;
        pTri_norms[i] = point3d{tri_norms[j], tri_norms[j+1], tri_norms[j+2]};
        pTri_barycs[i] = steps::math::tri_barycenter(pVerts[tri[0]], pVerts[tri[1]], pVerts[tri[2]]);
    }
    pTri_areas = tri_areas;
    for (auto a: pTri_areas) if (a<=0) throw ArgErr("triangle with non-positive area");
    pTri_diffboundaries.assign(pTrisN, nullptr);
    pTri_patches.assign(pTrisN, nullptr);

    // use tri data to make pBars, pTriBars.
    buildBarData();

    // copy tet information and compute tet barycentres.
    pTets.resize(pTetsN);
    pTet_barycentres.resize(pTetsN);
    for (uint i = 0, j = 0; i < pTetsN; ++i, j += 4) {
        tet_verts tet = {tets[j], tets[j+1], tets[j+2], tets[j+3]};

        pTets[i] = tet;
        pTet_barycentres[i] = steps::math::tet_barycenter(pVerts[tet[0]], pVerts[tet[1]], pVerts[tet[2]], pVerts[tet[3]]);
    }
    pTet_vols = tet_vols;
    for (auto v: pTet_vols) if (v<=0) throw ArgErr("tetrahedron with non-positive volume");

    // copy tet/tri and tet/tet neighbourhood infoirmation.
    pTri_tet_neighbours.resize(pTrisN);
    for (uint i = 0, j = 0; i < pTrisN; ++i, j += 2)
        pTri_tet_neighbours[i] = tri_tets{tri_tet_neighbs[j], tri_tet_neighbs[j+1]};

    pTet_tri_neighbours.resize(pTetsN);
    pTet_tet_neighbours.resize(pTetsN);
    for (uint i = 0, j = 0; i < pTetsN; ++i, j += 4) {
        pTet_tri_neighbours[i] = tet_tris{tet_tri_neighbs[j], tet_tri_neighbs[j+1], tet_tri_neighbs[j+2], tet_tri_neighbs[j+3]};
        pTet_tet_neighbours[i] = tet_tets{tet_tet_neighbs[j], tet_tet_neighbs[j+1], tet_tet_neighbs[j+2], tet_tet_neighbs[j+3]};
    }

    pTet_comps.assign(pTetsN, nullptr);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::Tetmesh::~Tetmesh(void)
{
    for (auto &membs: pMembs) delete membs.second;
    for (auto &diffb: pDiffBoundaries) delete diffb.second;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getVertex(uint vidx) const
{
    if (vidx >= pVertsN) throw steps::ArgErr("Vertex index is out of range.");
    return as_vector(pVerts[vidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getBar(uint bidx) const
{
    if (bidx >= pBarsN) throw steps::ArgErr("Bar index is out of range.");
    return as_vector(pBars[bidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getTri(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return as_vector(pTris[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTriArea(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return pTri_areas[tidx];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTriBarycenter(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return as_vector(pTri_barycs[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTriNorm(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return as_vector(pTri_norms[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch * stetmesh::Tetmesh::getTriPatch(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return pTri_patches[tidx];
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTriPatch(uint tidx, stetmesh::TmPatch * patch)
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    pTri_patches[tidx] = patch;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTriDiffBoundary(uint tidx, stetmesh::DiffBoundary * diffb)
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    pTri_diffboundaries[tidx] = diffb;
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::DiffBoundary * stetmesh::Tetmesh::getTriDiffBoundary(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return pTri_diffboundaries[tidx];
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setBarSDiffBoundary(uint bidx, stetmesh::SDiffBoundary * sdiffb)
{
    if (bidx >= pBarsN) throw steps::ArgErr("Bar index is out of range.");
    pBar_sdiffboundaries[bidx] = sdiffb;
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::SDiffBoundary * stetmesh::Tetmesh::getBarSDiffBoundary(uint bidx) const
{
    if (bidx >= pBarsN) throw steps::ArgErr("Bar index is out of range.");
    return pBar_sdiffboundaries[bidx];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getTriBars(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return as_vector(pTri_bars[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tetmesh::getTriTetNeighb(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");
    return as_vector(pTri_tet_neighbours[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tetmesh::getTriTriNeighb(uint tidx, const stetmesh::TmPatch * tmpatch) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");

    // Triangles are neighbours if they share a bar
    std::vector<int> neighbours(3,-1);
    tri_bars bars = pTri_bars[tidx];

    for (uint tri = 0; tri < pTrisN; ++tri) {
        if (tri == tidx || pTri_patches[tri] != tmpatch) continue;

        tri_bars neighbtribars = pTri_bars[tri];
        int next = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (neighbtribars[j] != bars[i]) continue;

                if (neighbours[i] != -1) {
                    std::ostringstream os;
                    os << "Error in Patch initialisation for '" << tmpatch->getID()
                       << "'. Patch triangle idx " << tidx << " found to have more than 3 neighbours.";
                    throw steps::ArgErr(os.str());
                }

                neighbours[i] = tri;
                next = 1;
                break;
            }
            if (next) break;
        }
    }
    return neighbours;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<int> stetmesh::Tetmesh::getTriTriNeighb(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");

    // Triangles are neighbours if they share a bar
    std::vector<int> neighbours(3,-1);
    tri_bars bars = pTri_bars[tidx];

    for (uint tri = 0; tri < pTrisN; ++tri) {
        if (tri == tidx || pTri_patches[tri] == nullptr) continue;

        tri_bars neighbtribars = pTri_bars[tri];
        int next = 0;
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (neighbtribars[j] != bars[i]) continue;

                if (neighbours[i] != -1) {
                    std::ostringstream os;
                    os << "Error: Triangle idx " << tidx << " found to have more than 3 neighbours.";
                    throw steps::ArgErr(os.str());
                }

                neighbours[i] = tri;
                next = 1;
                break;
            }
            if (next) break;
        }
    }
    return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

std::set<uint> stetmesh::Tetmesh::getTriTriNeighbs(uint tidx) const
{
    if (tidx >= pTrisN) throw steps::ArgErr("Triangle index is out of range.");

    // Triangles are neighbours if they share a bar

    std::set<uint> neighbours;
    tri_bars bars = pTri_bars[tidx];

    for (uint tri = 0; tri < pTrisN; ++tri) {
        if (tri == tidx) continue;

        tri_bars neighbtribars = pTri_bars[tri];
        for (int i = 0; i < 3; ++i) {
            for (int j = 0; j < 3; ++j) {
                if (neighbtribars[j] != bars[i]) continue;

                neighbours.insert(tri);
                goto next_tri;
            }
        }
next_tri: ;
    }

    return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

std::set<uint> stetmesh::Tetmesh::getBarTriNeighbs(uint bidx) const
{
    if (bidx >= pBarsN) throw steps::ArgErr("Bar index is out of range.");

    std::set<uint> neighbours;

    for (uint tri = 0; tri < pTrisN; ++tri) {
        tri_bars neighbtribars = pTri_bars[tri];
        for (int i = 0; i < 3; ++i) {
            if (neighbtribars[i] != bidx) continue;

            neighbours.insert(tri);
            goto next_tri;
        }
next_tri: ;
    }

    return neighbours;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setBarTris(uint bidx, int itriidx, int otriidx)
{
    if (bidx >= pBarsN) throw steps::ArgErr("Bar index is out of range.");
    if (itriidx < 0 or otriidx < 0 or itriidx >= pTrisN or otriidx >= pTrisN)
    {
    	throw steps::ArgErr("Invalid triangle index.");
   	}
    if (pBar_tri_neighbours[bidx][0] >= 0 or pBar_tri_neighbours[bidx][1] >= 0)
    {
    	std::ostringstream os;
    	os << "Bar " << bidx << " is part of more than one surface diffusion boundary.";
    	throw steps::ArgErr(os.str());
    }

    pBar_tri_neighbours[bidx][0] = itriidx;
    pBar_tri_neighbours[bidx][1] = otriidx;

}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_flipTriTetNeighb(uint tidx)
{
    assert(tidx < pTrisN);

    tri_tets &tt = pTri_tet_neighbours[tidx];
    std::swap(tt[0],tt[1]);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_flipTriVerts(uint tidx)
{
    assert(tidx < pTrisN);

    tri_verts &tri = pTris[tidx];
    std::swap(tri[0],tri[1]);
 
    // also recalculate the triangle's normal; the long way (though should be
    // negative of previous normal)

    pTri_norms[tidx] = steps::math::tri_normal(pVerts[tri[0]], pVerts[tri[1]], pVerts[tri[2]]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint>  stetmesh::Tetmesh::getTet(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return as_vector(pTets[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getMeshVolume() const
{
    return std::accumulate(pTet_vols.begin(),pTet_vols.end(),0.0);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTetVol(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return pTet_vols[tidx];
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getTetQualityRER(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");

    tet_verts verts = pTets[tidx];
    const point3d &v0 = pVerts[verts[0]], &v1 = pVerts[verts[1]],
                  &v2 = pVerts[verts[2]], &v3 = pVerts[verts[3]];

    return steps::math::tet_circumrad(v0, v1, v2, v3) / steps::math::tet_shortestedge(v0, v1, v2, v3);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getTetBarycenter(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return as_vector(pTet_barycentres[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmComp * stetmesh::Tetmesh::getTetComp(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return pTet_comps[tidx];
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::setTetComp(uint tidx, stetmesh::TmComp * comp)
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    pTet_comps[tidx] = comp;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getTetTriNeighb(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return as_vector(pTet_tri_neighbours[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<int> stetmesh::Tetmesh::getTetTetNeighb(uint tidx) const
{
    if (tidx >= pTetsN) throw steps::ArgErr("Tetrahedron index is out of range.");
    return as_vector(pTet_tet_neighbours[tidx]);
}

////////////////////////////////////////////////////////////////////////////////

int stetmesh::Tetmesh::findTetByPoint(std::vector<double> p) const
{
    point3d x{p[0],p[1],p[2]};
    if (!pBBox.contains(x)) return -1;

    for (int tidx = 0; tidx < pTetsN; ++tidx) {
        tet_verts v = pTets[tidx];
        if (steps::math::tet_inside(pVerts[v[0]], pVerts[v[1]], pVerts[v[2]], pVerts[v[3]], x))
            return tidx;
    }

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBoundMin(void) const
{
    return as_vector(pBBox.min());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBoundMax(void) const
{
    return as_vector(pBBox.max());
}

////////////////////////////////////////////////////////////////////////////////

// Created by weiliang 2010.02.02
std::vector<int> stetmesh::Tetmesh::getSurfTris(void) const
{
    std::vector<int> tribounds;
    for (int t = 0; t < pTrisN; t++) {
        tri_tets trineighbor = pTri_tet_neighbours[t];
        if (trineighbor[0] == -1 || trineighbor[1] == -1) {
            tribounds.push_back(t);
        }
    }
    return tribounds;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_checkMembID(std::string const & id) const
{
    checkID(id);
    if (pMembs.find(id) != pMembs.end())
        throw steps::ArgErr("'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleMembIDChange(std::string const & o, std::string const & n)
{
    MembPMapCI m_old = pMembs.find(o);
    assert(m_old != pMembs.end());

    if (o == n) return;
    _checkMembID(n);

    Memb * m = m_old->second;
    assert(m != 0);
    pMembs.erase(m->getID());                        // or s_old->first
    pMembs.insert(MembPMap::value_type(n, m));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleMembAdd(stetmesh::Memb * memb)
{
    assert(memb->getContainer() == this);
    _checkMembID(memb->getID());
    pMembs.insert(MembPMap::value_type(memb->getID(), memb));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleMembDel(stetmesh::Memb * memb)
{
    assert(memb->getContainer() == this);
    pMembs.erase(memb->getID());
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_checkDiffBoundaryID(std::string const & id) const
{
    checkID(id);
    if (pDiffBoundaries.find(id) != pDiffBoundaries.end())
        throw steps::ArgErr("'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleDiffBoundaryIDChange(std::string const & o, std::string const & n)
{
    DiffBoundaryPMapCI db_old = pDiffBoundaries.find(o);
    assert(db_old != pDiffBoundaries.end());

    if (o == n) return;
    _checkDiffBoundaryID(n);

    DiffBoundary * db = db_old->second;
    assert(db != 0);
    pDiffBoundaries.erase(db->getID());                        // or s_old->first
    pDiffBoundaries.insert(DiffBoundaryPMap::value_type(n, db));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleDiffBoundaryAdd(stetmesh::DiffBoundary * diffb)
{
    assert(diffb->getContainer() == this);
    _checkDiffBoundaryID(diffb->getID());
    pDiffBoundaries.insert(DiffBoundaryPMap::value_type(diffb->getID(), diffb));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleDiffBoundaryDel(stetmesh::DiffBoundary * diffb)
{
    assert(diffb->getContainer() == this);
    pDiffBoundaries.erase(diffb->getID());
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_checkSDiffBoundaryID(std::string const & id) const
{
    checkID(id);
    if (pSDiffBoundaries.find(id) != pSDiffBoundaries.end())
        throw steps::ArgErr("'" + id + "' is already in use.");
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleSDiffBoundaryIDChange(std::string const & o, std::string const & n)
{
    SDiffBoundaryPMapCI db_old = pSDiffBoundaries.find(o);
    assert(db_old != pSDiffBoundaries.end());

    if (o == n) return;
    _checkSDiffBoundaryID(n);

    SDiffBoundary * db = db_old->second;
    assert(db != 0);
    pSDiffBoundaries.erase(db->getID());                        // or s_old->first
    pSDiffBoundaries.insert(SDiffBoundaryPMap::value_type(n, db));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleSDiffBoundaryAdd(stetmesh::SDiffBoundary * sdiffb)
{
    assert(sdiffb->getContainer() == this);
    _checkSDiffBoundaryID(sdiffb->getID());
    pSDiffBoundaries.insert(SDiffBoundaryPMap::value_type(sdiffb->getID(), sdiffb));
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::_handleSDiffBoundaryDel(stetmesh::SDiffBoundary * sdiffb)
{
    assert(sdiffb->getContainer() == this);
    pSDiffBoundaries.erase(sdiffb->getID());
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::SDiffBoundary * stetmesh::Tetmesh::_getSDiffBoundary(uint gidx) const
{
    assert (gidx < pSDiffBoundaries.size());
    std::map<std::string, SDiffBoundary *>::const_iterator db_it = pSDiffBoundaries.begin();
    for (uint i=0; i< gidx; ++i) ++db_it;
    return db_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::Memb * stetmesh::Tetmesh::_getMemb(uint gidx) const
{
    assert (gidx < pMembs.size());
    std::map<std::string, Memb *>::const_iterator mb_it = pMembs.begin();
    for (uint i=0; i< gidx; ++i) ++mb_it;
    return mb_it->second;
}

////////////////////////////////////////////////////////////////////////////////

steps::tetmesh::DiffBoundary * stetmesh::Tetmesh::_getDiffBoundary(uint gidx) const
{
    assert (gidx < pDiffBoundaries.size());
    std::map<std::string, DiffBoundary *>::const_iterator db_it = pDiffBoundaries.begin();
    for (uint i=0; i< gidx; ++i) ++db_it;
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
void batch_copy_components_n(const std::vector<T> &items, I idx_iter, size_t n, J out_iter) {
    size_t index = 0;
    try {
        for (size_t i=0; i < n; ++i) {
            index = *idx_iter++;
            const auto &item = items.at(index);
            out_iter = std::copy(item.begin(),item.end(),out_iter);
        }
    }
    catch (std::out_of_range &e) {
        throw steps::ArgErr("Index out of range: no item with index "+std::to_string(index)+".");
    }
}

/* Copy n elements of the vector items to out_iter, where the elements
 * are specified by indices given by idx_iter.
 */

template <typename T, typename I, typename J>
void batch_copy_n(const std::vector<T> &items, I idx_iter, size_t n, J out_iter) {
    size_t index = 0;
    try {
        for (size_t i=0; i < n; ++i) {
            index = *idx_iter++;
            *out_iter++ = items.at(index);
        }
    }
    catch (std::out_of_range &e) {
        throw steps::ArgErr("Index out of range: no item with index "+std::to_string(index)+".");
    }
}
////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBatchTetBarycentres(std::vector<uint> const & tets) const
{
    uint ntets = tets.size();
    std::vector<double> data(ntets * 3);

    batch_copy_components_n(pTet_barycentres, tets.begin(), ntets, data.begin());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTetBarycentresNP(const unsigned int* indices, int input_size, double* centres, int output_size) const
{
    if (input_size * 3 != output_size)
        throw steps::ArgErr("Length of output array should be 3 * length of input array.");
    
    batch_copy_components_n(pTet_barycentres, indices, input_size, centres);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBatchTriBarycentres(std::vector<uint> const & tris) const
{
    uint ntris = tris.size();
    std::vector<double> data(ntris * 3);
    
    batch_copy_components_n(pTri_barycs, tris.begin(), ntris, data.begin());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTriBarycentresNP(const unsigned int* indices, int input_size, double* centres, int output_size) const
{
    if (input_size * 3 != output_size)
        throw steps::ArgErr("Length of output array should be 3 * length of input array.");
    
    batch_copy_components_n(pTri_barycs, indices, input_size, centres);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getBatchVertices(std::vector<uint> const & verts) const
{
    uint nverts = verts.size();
    std::vector<double> data(nverts * 3);
    
    batch_copy_components_n(pVerts, verts.begin(), nverts, data.begin());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchVerticesNP(const unsigned int* indices, int input_size, double* coordinates, int output_size) const
{
    if (input_size * 3 != output_size)
        throw steps::ArgErr("Length of output array (coordinates) should be 3 * length of input array (indices).");
    
    batch_copy_components_n(pVerts, indices, input_size, coordinates);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getBatchTris(std::vector<uint> const & tris) const
{
    uint ntris = tris.size();
    std::vector<uint> data(ntris * 3, 0);
    
    batch_copy_components_n(pTris, tris.begin(), ntris, data.begin());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTrisNP(const unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const
{
    if (input_size * 3 != output_size)
        throw steps::ArgErr("Length of output array should be 3 * length of input array.");
    
    batch_copy_components_n(pTris, t_indices, input_size, v_indices);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getBatchTets(std::vector<uint> const & tets) const
{
    uint ntets = tets.size();
    std::vector<uint> data(ntets * 4, 0);
    
    batch_copy_components_n(pTets, tets.begin(), ntets, data.begin());
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTetsNP(const unsigned int* t_indices, int input_size, unsigned int* v_indices, int output_size) const
{
    if (input_size * 4 != output_size)
        throw steps::ArgErr("Length of output array should be 4 * length of input array.");
    
    batch_copy_components_n(pTets, t_indices, input_size, v_indices);
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tetmesh::getTriVerticesSetSizeNP(const unsigned int* t_indices, int input_size) const
{
    std::set<uint> unique_indices;
    batch_copy_components_n(pTris, t_indices, input_size, std::inserter(unique_indices,unique_indices.begin()));

    return unique_indices.size();
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tetmesh::getTetVerticesSetSizeNP(const unsigned int* t_indices, int input_size) const
{
    std::set<uint> unique_indices;
    batch_copy_components_n(pTets, t_indices, input_size, std::inserter(unique_indices,unique_indices.begin()));
    
    return unique_indices.size();
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getTriVerticesMappingSetNP(const unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const
{
    if (input_size * 3 != t_vertices_size)
        throw steps::ArgErr("Length of t_vertices array should be 3 * length of input array.");
    
    // put unique vertices in v_set
    auto v_set_indices = make_unique_indexer<unsigned int>(v_set);
    int t_vertices_index = 0;

    for (int t = 0; t < input_size; ++t) {
        uint tidx = t_indices[t];
        if (tidx >= pTrisN) 
            throw steps::ArgErr("Index out of range: no item with index "+std::to_string(tidx)+".");
        
        for (unsigned int vidx: pTris[tidx]) 
            t_vertices[t_vertices_index++] = v_set_indices[vidx];
    }

    assert((int)v_set_indices.size() == v_set_size);
    assert(t_vertices_index == t_vertices_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getTetVerticesMappingSetNP(const unsigned int* t_indices, int input_size, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const
{
    if (input_size * 4 != t_vertices_size)
        throw steps::ArgErr("Length of t_vertices array should be 4 * length of input array.");
    
    auto v_set_indices = make_unique_indexer<unsigned int>(v_set);
    int t_vertices_index = 0;

    for (int t = 0; t < input_size; ++t) {
        uint tidx = t_indices[t];
        if (tidx >= pTetsN) 
            throw steps::ArgErr("Index out of range: no item with index "+std::to_string(tidx)+".");
        
        for (unsigned int vidx: pTets[tidx]) 
            t_vertices[t_vertices_index++] = v_set_indices[vidx];
    }
    
    assert((int)v_set_indices.size() == v_set_size);
    assert(t_vertices_index == t_vertices_size);
}

////////////////////////////////////////////////////////////////////////////////


void stetmesh::Tetmesh::genPointsInTet(unsigned tidx, unsigned npnts, double* coords, int coord_size) const
{
    if (npnts * 3 != coord_size)
        throw steps::ArgErr("Coordinate array size should be 3 * npnts.");

    if (tidx >= pTetsN) 
        throw steps::ArgErr("Index out of range: no tetrahedron with index "+std::to_string(tidx)+".");
    
    auto tet = pTets[tidx];
    point3d v[4] = {pVerts[tet[0]], pVerts[tet[1]], pVerts[tet[2]], pVerts[tet[3]]};
    
    for (uint p = 0, j = 0; p < npnts; p++, j+=3) {
        double s = ((double) rand() / (RAND_MAX));
        double t = ((double) rand() / (RAND_MAX));
        double u = ((double) rand() / (RAND_MAX));

        point3d x = steps::math::tet_ranpnt(v[0], v[1], v[2], v[3], s, t, u);
        coords[j]   = x[0];
        coords[j+1] = x[1];
        coords[j+2] = x[2];
    }
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::genPointsInTri(unsigned tidx, unsigned npnts, double* coords, int coord_size) const
{
    if (npnts * 3 != coord_size)
        throw steps::ArgErr("Coordinate array size should be 3 * npnts.");

    if (tidx >= pTrisN)
        throw steps::ArgErr("Index out of range: no triangle with index "+std::to_string(tidx)+".");
    
    auto tri = pTris[tidx];
    point3d v[3] = {pVerts[tri[0]], pVerts[tri[1]], pVerts[tri[2]]};
    
    for (uint p = 0, j = 0; p < npnts; p++, j+=3) {
        double s = ((double) rand() / (RAND_MAX));
        double t = ((double) rand() / (RAND_MAX));
        
        point3d x = steps::math::tri_ranpnt(v[0], v[1], v[2], s, t);
        coords[j]   = x[0];
        coords[j+1] = x[1];
        coords[j+2] = x[2];
    }
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::genTetVisualPointsNP(const unsigned int* indices, int index_size, const unsigned int* point_counts, int count_size, double* coords, int coord_size) const
{
    if (index_size != count_size)
        throw steps::ArgErr("Length of point_counts array should be length of indices array.");
    
    uint cpos = 0;
    for (uint i = 0; i < index_size; i++) {
        uint tidx = indices[i];
        uint npnts = point_counts[i];
        uint ncoords = 3*npnts;

        if (cpos + ncoords > coord_size) throw steps::ArgErr("Length of coords array too short.");

        genPointsInTet(indices[i], npnts, coords+cpos, ncoords);
        cpos += ncoords;
    }
    
    if (cpos != coord_size) throw steps::ArgErr("Length of coords array longer than expected.");
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::genTriVisualPointsNP(const unsigned int* indices, int index_size, const unsigned int* point_counts, int count_size, double* coords, int coord_size) const
{
    if (index_size != count_size)
        throw steps::ArgErr("Length of point_counts array should be length of indices array.");
    
    uint cpos = 0;
    for (uint i = 0; i < index_size; i++) {
        uint tidx = indices[i];
        uint npnts = point_counts[i];
        uint ncoords = 3*npnts;

        if (cpos + ncoords > coord_size) throw steps::ArgErr("Length of coords array too short.");

        genPointsInTri(indices[i], npnts, coords+cpos, ncoords);
        cpos += ncoords;
    }
    
    if (cpos != coord_size) throw steps::ArgErr("Length of coords array longer than expected.");
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTetVolsNP(const unsigned int* indices, int index_size, double* volumes, int volume_size) const
{
    if (index_size != volume_size)
        throw steps::ArgErr("Length of volumes array should be length of indices array.");
    
    batch_copy_n(pTet_vols, indices, index_size, volumes);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getBatchTriAreasNP(const unsigned int* indices, int index_size, double* areas, int area_size) const
{
    if (index_size != area_size)
        throw steps::ArgErr("Length of areas array should be length of indices array.");
    
    batch_copy_n(pTri_areas, indices, index_size, areas);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::reduceBatchTetPointCountsNP(const unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density)
{
    if (index_size != count_size)
        throw steps::ArgErr("Length of point_counts array should be length of indices array.");
    
    for (uint i = 0; i < index_size; i++) {
        uint tidx = indices[i];
        
        if (tidx >= pTetsN) 
            throw steps::ArgErr("Index out of range: no tetrahedron with index "+std::to_string(tidx)+".");
    
        point_counts[i] = std::min(point_counts[i], (uint)(max_density * pTet_vols[tidx]));
    }
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::reduceBatchTriPointCountsNP(const unsigned int* indices, int index_size, unsigned int* point_counts, int count_size, double max_density)
{
    if (index_size != count_size)
        throw steps::ArgErr("Length of point_counts array should be length of indices array.");
    
    for (uint i = 0; i < index_size; i++) {
        uint tidx = indices[i];
        
        if (tidx >= pTrisN) 
            throw steps::ArgErr("Index out of range: no triangle with index "+std::to_string(tidx)+".");
    
        point_counts[i] = std::min(point_counts[i], (uint)(max_density * pTri_areas[tidx]));
    }
}


////////////////////////////////////////////////////////////////////////
// ROI (Region of Interest) Data
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::addROI(std::string id, steps::tetmesh::ElementType type, std::set<uint> const &indices)
{
    if (mROI.find(id) != mROI.end()) {
        std::cerr << "Warning: ROI data with id " << id << " already exists. Use replaceROI() to replace the data.\n";
    }
    else
    {
        ROISet data(type, indices);
        mROI[id] = data;
    }
}

////////////////////////////////////////////////////////////////////////////////
void stetmesh::Tetmesh::removeROI(std::string id)
{
    std::map<std::string, ROISet>::iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
    }
    else
    {
        mROI.erase(it);
    }
}

////////////////////////////////////////////////////////////////////////////////
void stetmesh::Tetmesh::replaceROI(std::string id, steps::tetmesh::ElementType type, std::set<uint> const &indices)
{
    std::map<std::string, ROISet>::iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
    }
    else
    {
        it->second.type = type;
        it->second.indices.assign(indices.begin(), indices.end());
    }
}

////////////////////////////////////////////////////////////////////////////////
steps::tetmesh::ElementType stetmesh::Tetmesh::getROIType(std::string id) const
{
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return steps::tetmesh::ELEM_UNDEFINED;
    }

    return it->second.type;
}

////////////////////////////////////////////////////////////////////////////////
std::vector<uint> stetmesh::Tetmesh::getROIData(std::string id) const
{
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return std::vector<uint> ();
    }
    return it->second.indices;
}

////////////////////////////////////////////////////////////////////////////////
uint stetmesh::Tetmesh::getROIDataSize(std::string id) const
{
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return ELEM_UNDEFINED;
    }

    return it->second.indices.size();
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tetmesh::getNROIs(void)
{
    return mROI.size();
}

////////////////////////////////////////////////////////////////////////////////

stetmesh::ROISet stetmesh::Tetmesh::getROI(std::string id) const
{
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return stetmesh::ROISet();
    }

    return it->second;
}


////////////////////////////////////////////////////////////////////////////////

std::vector<std::string> stetmesh::Tetmesh::getAllROINames(void)
{
    std::vector<std::string> output;
    
    for (std::map<std::string, ROISet>::iterator it = mROI.begin(); it != mROI.end(); it++) {
        output.push_back(it->first);
    }
    return output;
}

////////////////////////////////////////////////////////////////////////////////

uint* stetmesh::Tetmesh::_getROIData(std::string id) const
{
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return NULL;
    }

    return const_cast<uint*>(&(it->second.indices.front()));
}

////////////////////////////////////////////////////////////////////////////////

bool stetmesh::Tetmesh::checkROI(std::string id, steps::tetmesh::ElementType type, uint count, bool warning) const
{
    bool condition = true;
    std::map<std::string, ROISet>::const_iterator it = mROI.find(id);
    if (warning && it == mROI.end()) {
        std::cerr << "Warning: Unable to find ROI data with id " << id << ".\n";
        return false;
    }
    if (warning && it->second.type != type) {
        std::cerr << "Warning: Element type mismatch for ROI " << id << ".\n";
        return false;
    }
    if (warning && count != 0 && it->second.indices.size() != count) {
        std::cerr << "Warning: Element count mismatch for ROI " << id << ".\n";
        return false;
    }
    return true;
}

////////////////////////////////////////////////////////////////////////
// ROI Data Access
////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getROITetBarycentres(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    std::vector<double> data(inputsize * 3, 0.0);
    
    getBatchTetBarycentresNP(indices, inputsize, const_cast<double*>(&data.front()), data.size());

    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITetBarycentresNP(std::string ROI_id, double* centres, int output_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, output_size / 3)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTetBarycentresNP(indices, inputsize, centres, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getROITriBarycentres(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    std::vector<double> data(inputsize * 3, 0.0);
    
    getBatchTriBarycentresNP(indices, inputsize, const_cast<double*>(&data.front()), data.size());
    
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITriBarycentresNP(std::string ROI_id, double* centres, int output_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, output_size / 3)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTriBarycentresNP(indices, inputsize, centres, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::Tetmesh::getROIVertices(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_VERTEX)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    std::vector<double> data(inputsize * 3, 0.0);
    getBatchVerticesNP(indices, inputsize, const_cast<double*>(&data.front()), data.size());
    
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROIVerticesNP(std::string ROI_id, double* coordinates, int output_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_VERTEX, output_size / 3)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchVerticesNP(indices, inputsize, coordinates, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getROITris(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    std::vector<uint> data(inputsize * 3, 0.0);
    getBatchTrisNP(indices, inputsize, const_cast<uint*>(&data.front()), data.size());
    
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITrisNP(std::string ROI_id, unsigned int* v_indices, int output_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, output_size / 3)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTrisNP(indices, inputsize, v_indices, output_size);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> stetmesh::Tetmesh::getROITets(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    std::vector<uint> data(inputsize * 4, 0.0);
    getBatchTetsNP(indices, inputsize, const_cast<uint*>(&data.front()), data.size());
    
    return data;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITetsNP(std::string ROI_id, unsigned int* v_indices, int output_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, output_size / 4)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTetsNP(indices, inputsize, v_indices, output_size);
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tetmesh::getROITriVerticesSetSizeNP(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    return getTriVerticesSetSizeNP(indices, inputsize);
}

////////////////////////////////////////////////////////////////////////////////

uint stetmesh::Tetmesh::getROITetVerticesSetSizeNP(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    return getTetVerticesSetSizeNP(indices, inputsize);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITriVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, t_vertices_size / 3)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getTriVerticesMappingSetNP(indices, inputsize, t_vertices, t_vertices_size, v_set, v_set_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITetVerticesMappingSetNP(std::string ROI_id, unsigned int* t_vertices, int t_vertices_size, unsigned int* v_set, int v_set_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, t_vertices_size / 4)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getTriVerticesMappingSetNP(indices, inputsize, t_vertices, t_vertices_size, v_set, v_set_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::genROITetVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, count_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    genTetVisualPointsNP(indices, inputsize, point_counts, count_size, coords, coord_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::genROITriVisualPointsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double* coords, int coord_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, count_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    genTriVisualPointsNP(indices, inputsize, point_counts, count_size, coords, coord_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITetVolsNP(std::string ROI_id, double* volumes, int volume_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, volume_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTetVolsNP(indices, inputsize, volumes, volume_size);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::getROITriAreasNP(std::string ROI_id, double* areas, int area_size) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, area_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    getBatchTriAreasNP(indices, inputsize, areas, area_size);
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getROIVol(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET)) throw steps::ArgErr();
    uint *indices = _getROIData(ROI_id);
    int roisize = getROIDataSize(ROI_id);
    
    double sum_vol = 0.0;
    for (uint t = 0; t < roisize; t++) {
        sum_vol += getTetVol(indices[t]);
    }
    
    return sum_vol;
}

////////////////////////////////////////////////////////////////////////////////

double stetmesh::Tetmesh::getROIArea(std::string ROI_id) const
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI)) throw steps::ArgErr();
    uint *indices = _getROIData(ROI_id);
    int roisize = getROIDataSize(ROI_id);
    
    double sum_area = 0.0;
    for (uint t = 0; t < roisize; t++) {
        sum_area += getTriArea(indices[t]);
    }
    
    return sum_area;
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::reduceROITetPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density)
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TET, count_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    reduceBatchTetPointCountsNP(indices, inputsize, point_counts, count_size, max_density);
}

////////////////////////////////////////////////////////////////////////////////

void stetmesh::Tetmesh::reduceROITriPointCountsNP(std::string ROI_id, unsigned int* point_counts, int count_size, double max_density)
{
    if (!checkROI(ROI_id, steps::tetmesh::ELEM_TRI, count_size)) throw steps::ArgErr();
    
    uint *indices = _getROIData(ROI_id);
    int inputsize = getROIDataSize(ROI_id);
    
    reduceBatchTriPointCountsNP(indices, inputsize, point_counts, count_size, max_density);
}

// END
