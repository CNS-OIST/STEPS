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

// Standard headers.
#include <cassert>
#include <algorithm>
#include <unordered_set>
#include <vector>
#include <sstream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/math/point.hpp"
#include "steps/util/collections.hpp"

namespace stetmesh = steps::tetmesh;

using steps::math::point3d;

////////////////////////////////////////////////////////////////////////////////

stetmesh::TmPatch::TmPatch(std::string const & id, Tetmesh * container,
                           std::vector<uint> const & tris, steps::wm::Comp* wmicomp,
                           steps::wm::Comp* wmocomp)
: steps::wm::Patch(id, container, wmicomp, wmocomp, 0.0)
, pTetmesh(container)
, pTrisN(0)
{
    if (pTetmesh == 0)
        throw steps::ArgErr("No mesh provided to Patch initializer function.");

    // upcast the compartment pointers for this overloaded constructor
    stetmesh::TmComp *icomp = dynamic_cast<stetmesh::TmComp *>(wmicomp);
    stetmesh::TmComp *ocomp = dynamic_cast<stetmesh::TmComp *>(wmocomp);

    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countTris() -1);

    // The patch's area - contributed to from all triangles
    double area = 0.0;

    std::unordered_set<uint> visited_tris(tris.size());
    for (uint tri: tris) {
        if (visited_tris.count(tri)) continue;
        visited_tris.insert(tri);

        if (tri > maxidx)
            throw steps::ArgErr("Invalid triangle index "+std::to_string(tri)+".");

        if (pTetmesh->getTriPatch(tri) != nullptr)
            throw steps::ArgErr("Triangle with index "+std::to_string(tri)+" already belongs to a patch.");

        if (pTetmesh->getTriDiffBoundary(tri) != nullptr)
            throw steps::ArgErr("Triangle with index "+std::to_string(tri)+" belongs to a diffusion boundary.");

        // Add triangle if compartments match those of patch, flipping
        // triangle neighbours if required.

        const int *tri_tets = pTetmesh->_getTriTetNeighb(tri);

        const stetmesh::TmComp *tri_comps[2] = {
            tri_tets[0]==-1? nullptr: pTetmesh->getTetComp(tri_tets[0]),
            tri_tets[1]==-1? nullptr: pTetmesh->getTetComp(tri_tets[1]),
        };

        if (tri_comps[1] == icomp && tri_comps[0] == ocomp) {
            pTetmesh->_flipTriTetNeighb(tri);
            pTri_indices.push_back(tri);
        }
        else if (tri_comps[0] == icomp && tri_comps[1] == ocomp) {
            pTri_indices.push_back(tri);
        }
        else throw steps::ArgErr("Triangle with index "+std::to_string(tri)
                +" has incompatible compartments for patch.");

        // If triangle normal does not point towards inner (first neighbour)
        // tetrahedron, flip triangle.

        // (sgy: Is it safe to assume inner tet is not -1? i.e. icomp never null?)
        assert(tri_tets[0] >= 0);

        point3d b_to_b = pTetmesh->_getTriBarycenter(tri) - pTetmesh->_getTetBarycenter(tri_tets[0]);
        if (dot(b_to_b, pTetmesh->_getTriNorm(tri)) < 0)
            pTetmesh->_flipTriVerts(tri);

        // Update area, patch, bounding box.
        area += pTetmesh->getTriArea(tri);
        pTetmesh->setTriPatch(tri, this);

        const uint *tri_verts = pTetmesh->_getTri(tri);
        for (uint j = 0; j < 3; ++j) {
            pBBox.insert(pTetmesh->_getVertex(tri_verts[j]));
        }
    } // end of loop over all tris (argument to constructor)

    pTrisN = pTri_indices.size();
    setArea(area);
}

std::vector<bool> stetmesh::TmPatch::isTriInside(const std::vector<uint> &tris) const
{
    return steps::util::map_membership(tris, pTri_indices);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmPatch::getBoundMin(void) const
{
    return steps::util::as_vector(pBBox.min());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> stetmesh::TmPatch::getBoundMax(void) const
{
    return steps::util::as_vector(pBBox.max());
}

////////////////////////////////////////////////////////////////////////////////

// END
