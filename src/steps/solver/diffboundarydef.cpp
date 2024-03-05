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

#include "diffboundarydef.hpp"

#include "geom/diffboundary.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

////////////////////////////////////////////////////////////////////////////////

DiffBoundarydef::DiffBoundarydef(Statedef&, diffboundary_global_id idx, tetmesh::DiffBoundary& db)
    : pIdx(idx)
    , pName(db.getID())
    , pTris(db._getAllTriIndices()) {
    const auto& comps = db.getComps();
    pCompA_temp = comps[0];
    pCompB_temp = comps[1];
    AssertLog(pCompA_temp != nullptr);
    AssertLog(pCompB_temp != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::checkpoint(std::fstream& /*cp_file*/) const {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::setup(Statedef& sd) {
    AssertLog(pSetupdone == false);

    pCompA = sd.getCompIdx(*pCompA_temp);
    pCompB = sd.getCompIdx(*pCompB_temp);
    pSetupdone = true;
}

}  // namespace steps::solver
