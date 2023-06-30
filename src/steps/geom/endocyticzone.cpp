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

// Standard headers.
#include <algorithm>
#include <string>
#include <vector>

// STEPS headers.
#include "geom/endocyticzone.hpp"
#include "geom/tmpatch.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
// logging
#include "easylogging++.h"

namespace steps::tetmesh {

EndocyticZone::EndocyticZone(std::string const& id,
                             TmPatch* patch,
                             std::vector<index_t> const& tris)
    : pID(id)
    , pPatch(patch) {
    ArgErrLogIf(pPatch == nullptr, "No patch provided to EndocyticZone initializer function.");

    ArgErrLogIf(tris.empty(), "The triangle list is empty.");

    const auto& trisInside = pPatch->isTriInside(tris);
    bool allTrisInPatch =
        std::all_of(trisInside.begin(), trisInside.end(), [](bool b) { return b; });

    ArgErrLogIf(not allTrisInPatch, "Some triangles are not part of the provided patch.");

    for (auto tri: tris) {
        pTri_indices.emplace_back(tri);
    }

    pPatch->_addEndocyticZone(this);
}

////////////////////////////////////////////////////////////////////////////////

EndocyticZone::~EndocyticZone() = default;

}  // namespace steps::tetmesh
