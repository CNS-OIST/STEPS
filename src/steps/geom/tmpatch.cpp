/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/math/point.hpp"
#include "steps/util/collections.hpp"
// logging
#include "easylogging++.h"

namespace steps {
namespace tetmesh {

using math::point3d;

////////////////////////////////////////////////////////////////////////////////

TmPatch::TmPatch(std::string const & id, Tetmesh * container,
                           std::vector<index_t> const & tris,
                           steps::wm::Comp* wmicomp,
                           steps::wm::Comp* wmocomp)
: steps::wm::Patch(id, container, wmicomp, wmocomp, 0.0)
, pTetmesh(container)
, pTrisN(0)
{
    if (pTetmesh == nullptr) {
        ArgErrLog("No mesh provided to Patch initializer function.");
    }

    // upcast the compartment pointers for this overloaded constructor

    // Note that for well-mixed compartment this will return nullptr
    auto icomp = dynamic_cast<TmComp *>(wmicomp);
    auto ocomp = dynamic_cast<TmComp *>(wmocomp);

    // The maximum triangle index in tetrahedral mesh
    auto maxidx = pTetmesh->countTris() - 1;

    // The patch's area - contributed to from all triangles
    double area = 0.0;

    std::unordered_set<index_t> visited_tris(tris.size());
    for (auto tri: tris) {
        if (visited_tris.count(tri) > 0) {
            continue;
        }
        visited_tris.insert(tri);

        if (tri > maxidx) {
            ArgErrLog("Invalid triangle index " + std::to_string(tri) + ".");
        }

        if (pTetmesh->getTriPatch(tri) != nullptr) {
            ArgErrLog("Triangle with index " + std::to_string(tri) + " already belongs to a patch.");
        }

        if (pTetmesh->getTriDiffBoundary(tri) != nullptr) {
            ArgErrLog("Triangle with index " + std::to_string(tri) + " belongs to a diffusion boundary.");
        }

        // Add triangle if compartments match those of patch, flipping
        // triangle neighbours if required.

        const auto *tri_tets = pTetmesh->_getTriTetNeighb(tri);

        // Weiliang: bug fix for trianglar patch with pure well-mixed compartments

        // this implementation guarantee the following:
        // 1. if a tri's tet neighbor is in inner comp, its id is in tri_tet[0]
        // 2. if a tri's tet neighbor is in outer comp, its id is in tri_tet[1]

        // note that if the tet neighbor is not in either comp (well-mixed compartment),
        // there is no guarantee that tri_tets[0] != -1

        // if tri_tet[0] doesn't exist but tri_tet[1] exists
        if (tri_tets[0] == UNKNOWN_TET && tri_tets[1] != UNKNOWN_TET) {
            auto tri_comp = pTetmesh->getTetComp(tri_tets[1]);
            // if tri_comp is inner
            // we flip it so that tri_tets[0] is now in inner comp
            if (tri_comp != nullptr && tri_comp == icomp) {
                pTetmesh->_flipTriTetNeighb(tri);
            }
        }
        else if (tri_tets[0] != UNKNOWN_TET && tri_tets[1] == UNKNOWN_TET) {
            auto tri_comp = pTetmesh->getTetComp(tri_tets[0]);
            // if it is outcomp or nullptr
            // we flip it so that now tri_tets[1] is outer comp
            if (tri_comp != nullptr && tri_comp == ocomp) {
                pTetmesh->_flipTriTetNeighb(tri);
            }
        }
        // both sides have tetrahedrons
        else if (tri_tets[0] != UNKNOWN_TET && tri_tets[1] != UNKNOWN_TET) {
            auto tri_comp0 = pTetmesh->getTetComp(tri_tets[0]);
            auto tri_comp1 = pTetmesh->getTetComp(tri_tets[1]);
            // if tri_comp0 is not nullptr (no well-mixed) and is outcomp,
            // we flip so that tri_tets[1] is now ocomp
            if (tri_comp0 != nullptr && tri_comp0 == ocomp) {
                pTetmesh->_flipTriTetNeighb(tri);
            }
            // if tri_comp1 is not nullptr (no well-mixed) and is innercomp,
            // we flip so that tri_tets[0] is now icomp
            else if (tri_comp1 != nullptr &&  tri_comp1 == icomp) {
                pTetmesh->_flipTriTetNeighb(tri);
            }
            // in other cases (including both comps are well-mixed), we don't need to flip

        }
        // both sides have no tet, this shouldn't happen
        else {
            ArgErrLog("Triangle with index "+std::to_string(tri)
                                +" has no neighboring tetrahedron.");
        }

        pTri_indices.push_back(tri);

        // from above we are sure that if updated_tri_tets are used in the simulation,
        // updated_tri_tets[0] will be inner tet, updated_tri_tets[1] will be outer tet.
        // note that they may not be used in the simulation, in which case flipping or not doesn't matter

        const auto *updated_tri_tets = pTetmesh->_getTriTetNeighb(tri);
        
        // tri_tets[0] exists
        if (updated_tri_tets[0] != UNKNOWN_TET) {
            point3d b_to_b = pTetmesh->_getTriBarycenter(tri) - pTetmesh->_getTetBarycenter(updated_tri_tets[0]);
            if (b_to_b.dot(pTetmesh->_getTriNorm(tri)) < 0) {
                pTetmesh->_flipTriVerts(tri);
            }
        }
        // tri_tets[1] exists
        else if (updated_tri_tets[1] != UNKNOWN_TET) {
            const point3d b_to_b = pTetmesh->_getTriBarycenter(tri) - pTetmesh->_getTetBarycenter(updated_tri_tets[1]);
            if (b_to_b.dot(pTetmesh->_getTriNorm(tri)) >= 0) {
                pTetmesh->_flipTriVerts(tri);
            }
        }

        // Update area, patch, bounding box.
        area += pTetmesh->getTriArea(tri);
        pTetmesh->setTriPatch(tri, this);

        const auto tri_verts = pTetmesh->_getTri(tri);
        for (uint j = 0; j < 3; ++j) {
            pBBox.insert(pTetmesh->_getVertex(tri_verts[j]));
        }
    } // end of loop over all tris (argument to constructor)

    pTrisN = pTri_indices.size();
    setArea(area);
}

std::vector<bool> TmPatch::isTriInside(const std::vector<index_t> &tris) const
{
    return steps::util::map_membership(tris, pTri_indices);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmPatch::getBoundMin() const
{
    return steps::util::as_vector(pBBox.min());
}

////////////////////////////////////////////////////////////////////////////////

std::vector<double> TmPatch::getBoundMax() const
{
    return steps::util::as_vector(pBBox.max());
}

} // namespace tetmesh
} // namespace steps
