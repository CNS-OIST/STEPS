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
#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/error.hpp"
#include "steps/geom/diffboundary.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/util/collections.hpp"

// logging
#include "easylogging++.h"

namespace steps {
namespace tetmesh {

DiffBoundary::DiffBoundary(std::string id, Tetmesh * container,
            std::vector<index_t> const & tris)
: pID(std::move(id))
, pTetmesh(container)
{
    if (pTetmesh == nullptr) {
        ArgErrLog("No mesh provided to Diffusion Boundary initializer function.");
    }

    if (tris.empty()) {
        ArgErrLog("No triangles provided to Diffusion Boundary initializer function.");
    }

    // The maximum triangle index in tetrahedral mesh
    const auto maxidx = (pTetmesh->countTris() -1);

    std::unordered_set<index_t> visited_tris(tris.size());
    for (auto tri: tris) {
        if (visited_tris.find(tri) != visited_tris.end()) {
            continue;
        }
        visited_tris.insert(tri);

        if (tri > maxidx) {
            ArgErrLog("Invalid triangle index " + std::to_string(tri) + ".");
        }

        if (pTetmesh->getTriPatch(tri) != nullptr) {
            ArgErrLog("Triangle with index " + std::to_string(tri) + " belongs to a patch.");
        }

        if (pTetmesh->getTriDiffBoundary(tri) != nullptr) {
            ArgErrLog("Triangle with index " + std::to_string(tri) + " already belongs to a diffusion boundary.");
        }

        const auto *tri_tets = pTetmesh->_getTriTetNeighb(tri);

        if (tri_tets[0] == UNKNOWN_TET || tri_tets[1] == UNKNOWN_TET) {
            ArgErrLog("Triangle with index " + std::to_string(tri) + " is on mesh surface.");
        }

        auto *tri_icomp = pTetmesh->getTetComp(tri_tets[0]);
        auto *tri_ocomp = pTetmesh->getTetComp(tri_tets[1]);

        if (tri_icomp == nullptr || tri_ocomp == nullptr) {
            ArgErrLog(
              "Triangle with index " + std::to_string(tri) + " does not have both an inner and outer compartment.");
        }

        if (pIComp == nullptr) {
            // Set diffboundary compartments from first tet.
            pIComp = tri_icomp;
            pOComp = tri_ocomp;
        }
        else {
            // Tet compartments may be other way around.
            if (pIComp != tri_icomp) {
                std::swap(tri_icomp,tri_ocomp);
            }

            if (pIComp != tri_icomp || pOComp != tri_ocomp) {
                ArgErrLog("Triangle with index " + std::to_string(tri) + " has incompatible compartments.");
            }
        }

        pTri_indices.push_back(tri);
        pTetmesh->setTriDiffBoundary(tri, this);
    } // end of loop over all tris (argument to constructor)

    pTrisN = pTri_indices.size();
    pTetmesh->_handleDiffBoundaryAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::setID(std::string const & id)
{
    AssertLog(pTetmesh != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pTetmesh->_handleDiffBoundaryIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> DiffBoundary::isTriInside(std::vector<index_t> const &tris) const
{
    return steps::util::map_membership(tris, pTri_indices);
}

} // namespace tetmesh
} // namespace steps
