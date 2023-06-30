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

// STEPS headers.
#include "diffboundarydef.hpp"
#include "geom/comp.hpp"
#include "types.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {


////////////////////////////////////////////////////////////////////////////////

DiffBoundarydef::DiffBoundarydef(Statedef* sd,
                                 diffboundary_global_id idx,
                                 tetmesh::DiffBoundary* db)
    : pStatedef(sd)
    , pIdx(idx)
    , pCompA_temp(nullptr)
    , pCompB_temp(nullptr) {
    AssertLog(pStatedef != nullptr);
    AssertLog(db != nullptr);

    pName = db->getID();
    pTris = db->_getAllTriIndices();
    std::vector<wm::Comp*> comps = db->getComps();
    pCompA_temp = comps[0];
    pCompB_temp = comps[1];
    AssertLog(pCompA_temp != nullptr);
    AssertLog(pCompB_temp != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

DiffBoundarydef::~DiffBoundarydef() = default;

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundarydef::setup() {
    AssertLog(pSetupdone == false);

    pCompA = pStatedef->getCompIdx(pCompA_temp);
    pCompB = pStatedef->getCompIdx(pCompB_temp);
    pSetupdone = true;
}

}  // namespace steps::solver
