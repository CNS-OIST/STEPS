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

#include "sdiffboundary.hpp"

#include <algorithm>
#include <cassert>
#include <sstream>
#include <string>
#include <unordered_set>
#include <vector>

#include <easylogging++.h>

#include "tmpatch.hpp"

#include "util/collections.hpp"
#include "util/error.hpp"

namespace steps::tetmesh {

////////////////////////////////////////////////////////////////////////////////

SDiffBoundary::SDiffBoundary(std::string id,
                             Tetmesh* container,
                             std::vector<index_t> const& bars,
                             std::vector<TmPatch*> const& patches)
    : pID(std::move(id))
    , pTetmesh(container) {
    ArgErrLogIf(pTetmesh == nullptr,
                "No mesh provided to Surface Diffusion Boundary initializer function.");

    ArgErrLogIf(bars.empty(),
                "No triangles provided to Surface Diffusion Boundary "
                "initializer function.");

    ArgErrLogIf(patches.size() != 2,
                "The number of patches provided to Surface Diffusion initializer "
                "function must be length 2.");

    // The maximum triangle index in tetrahedral mesh
    uint maxidx = (pTetmesh->countBars() - 1);

    std::unordered_set<bar_id_t> visited_bars(bars.size());
    for (uint bar_idx: bars) {
        bar_id_t bar(bar_idx);
        if (visited_bars.count(bar) != 0) {
            continue;
        }
        visited_bars.insert(bar);

        ArgErrLogIf(bar > maxidx, "Invalid bar index " + std::to_string(bar) + ".");

        ArgErrLogIf(pTetmesh->getBarSDiffBoundary(bar) != nullptr,
                    "Bar with index " + std::to_string(bar) +
                        " already belongs to a surface diffusion boundary.");

        // Need to find one triangle from inner patch and one from outer in this
        // lot:
        const std::set<triangle_global_id> bartris = pTetmesh->getBarTriNeighbs(bar);

        // Must be one and only one copy of inner patch and outer patch
        triangle_global_id innertriidx;
        triangle_global_id outertriidx;

        for (const auto& tri: bartris) {
            TmPatch* tri_patch = pTetmesh->getTriPatch(tri);

            if (tri_patch == patches[0]) {
                if (innertriidx.unknown()) {
                    innertriidx = tri;
                } else {
                    ArgErrLog("Duplicate copy of patch" + patches[0]->getID() +
                              " in connected triangles to bar " + std::to_string(bar));
                }
            } else if (tri_patch == patches[1]) {
                if (outertriidx.unknown()) {
                    outertriidx = tri;
                } else {
                    ArgErrLog("Duplicate copy of patch" + patches[1]->getID() +
                              " in connected triangles to bar " + std::to_string(bar));
                }
            }
        }

        ArgErrLogIf(innertriidx.unknown(),
                    "Patch" + patches[0]->getID() + " is not connected to bar " +
                        std::to_string(bar));

        ArgErrLogIf(outertriidx.unknown(),
                    "Patch" + patches[1]->getID() + " is not connected to bar " +
                        std::to_string(bar));

        pBar_indices.push_back(bar.get());
        pTetmesh->setBarSDiffBoundary(bar, this);

        pTetmesh->setBarTris(bar, innertriidx, outertriidx);

    }  // end of loop over all tris (argument to constructor)

    // Tests passed, last thing to do is record these patches
    pIPatch = patches[0];
    pOPatch = patches[1];

    pTetmesh->_handleSDiffBoundaryAdd(this);
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundary::setID(std::string const& id) {
    AssertLog(pTetmesh != nullptr);
    if (id == pID) {
        return;
    }
    // The following might raise an exception, e.g. if the new ID is not
    // valid or not unique. If this happens, we don't catch but simply let
    // it pass by into the Python layer.
    pTetmesh->_handleSDiffBoundaryIDChange(pID, id);
    // This line will only be executed if the previous call didn't raise
    // an exception.
    pID = id;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<bool> SDiffBoundary::isBarInside(const std::vector<index_t>& bars) const {
    return util::map_membership(bars, pBar_indices);
}

}  // namespace steps::tetmesh
