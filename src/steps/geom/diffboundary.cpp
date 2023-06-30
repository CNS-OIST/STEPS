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

#include "diffboundary.hpp"

#include <algorithm>
#include <string>
#include <unordered_set>
#include <vector>

#include <easylogging++.h>

#include "util/collections.hpp"
#include "util/error.hpp"

namespace steps::tetmesh {

DiffBoundary::DiffBoundary(std::string id, Tetmesh* container, std::vector<index_t> const& tris)
    : pID(std::move(id))
    , pTetmesh(container) {
    ArgErrLogIf(pTetmesh == nullptr,
                "No mesh provided to Diffusion Boundary initializer function.");
    ArgErrLogIf(tris.empty(), "No triangles provided to Diffusion Boundary initializer function.");

    // The maximum triangle index in tetrahedral mesh
    const auto maxidx = (pTetmesh->countTris() - 1);

    std::unordered_set<triangle_global_id> visited_tris(tris.size());
    for (auto tri_idx: tris) {
        triangle_global_id tri(tri_idx);
        if (visited_tris.find(tri) != visited_tris.end()) {
            continue;
        }
        visited_tris.insert(tri);

        ArgErrLogIf(tri > maxidx, "Invalid triangle index " + std::to_string(tri) + ".");
        ArgErrLogIf(pTetmesh->getTriPatch(tri) != nullptr,
                    "Triangle with index " + std::to_string(tri) + " belongs to a patch.");
        ArgErrLogIf(pTetmesh->getTriDiffBoundary(tri) != nullptr,
                    "Triangle with index " + std::to_string(tri) +
                        " already belongs to a diffusion boundary.");
        auto tri_tets = pTetmesh->_getTriTetNeighb(tri);
        ArgErrLogIf(tri_tets[0].unknown() || tri_tets[1].unknown(),
                    "Triangle with index " + std::to_string(tri) + " is on mesh surface.");

        auto* tri_icomp = pTetmesh->getTetComp(tri_tets[0]);
        auto* tri_ocomp = pTetmesh->getTetComp(tri_tets[1]);

        ArgErrLogIf(tri_icomp == nullptr || tri_ocomp == nullptr,
                    "Triangle with index " + std::to_string(tri) +
                        " does not have both an inner and outer compartment.");

        if (pIComp == nullptr) {
            // Set diffboundary compartments from first tet.
            pIComp = tri_icomp;
            pOComp = tri_ocomp;
        } else {
            // Tet compartments may be other way around.
            if (pIComp != tri_icomp) {
                std::swap(tri_icomp, tri_ocomp);
            }

            ArgErrLogIf(pIComp != tri_icomp || pOComp != tri_ocomp,
                        "Triangle with index " + std::to_string(tri) +
                            " has incompatible compartments.");
        }

        pTri_indices.push_back(tri);
        pTetmesh->setTriDiffBoundary(tri, this);
    }  // end of loop over all tris (argument to constructor)

    pTrisN = pTri_indices.size();
    pTetmesh->_handleDiffBoundaryAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::setID(std::string const& id) {
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

std::vector<bool> DiffBoundary::isTriInside(std::vector<index_t> const& tris) const {
    return util::map_membership(tris, pTri_indices);
}

}  // namespace steps::tetmesh
