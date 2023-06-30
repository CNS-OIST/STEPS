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
#include "sdiffboundarydef.hpp"
#include "geom/patch.hpp"
#include "types.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

SDiffBoundarydef::SDiffBoundarydef(Statedef* sd,
                                   sdiffboundary_global_id idx,
                                   tetmesh::SDiffBoundary* sdb)
    : pStatedef(sd)
    , pSetupdone(false)
    , pIdx(idx)
    , pPatchA_temp(nullptr)
    , pPatchB_temp(nullptr) {
    AssertLog(pStatedef != nullptr);
    AssertLog(sdb != nullptr);

    pName = sdb->getID();
    pBars = sdb->_getAllBarIndices();
    const auto& patches = sdb->getPatches();
    pPatchA_temp = patches[0];
    pPatchB_temp = patches[1];
    AssertLog(pPatchA_temp != nullptr);
    AssertLog(pPatchB_temp != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

SDiffBoundarydef::~SDiffBoundarydef() = default;

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundarydef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundarydef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundarydef::setup() {
    AssertLog(pSetupdone == false);

    pPatchA = pStatedef->getPatchIdx(pPatchA_temp);
    pPatchB = pStatedef->getPatchIdx(pPatchB_temp);
    pSetupdone = true;
}

}  // namespace steps::solver
