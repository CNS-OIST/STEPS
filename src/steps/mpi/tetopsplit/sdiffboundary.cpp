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
#include "sdiffboundary.hpp"
// logging
#include "util/error.hpp"
#include <easylogging++.h>

namespace steps::mpi::tetopsplit {

SDiffBoundary::SDiffBoundary(solver::SDiffBoundarydef* sdbdef)
    : pSDiffBoundarydef(sdbdef)
    , pSetPatches(false)
    , pPatchA(nullptr)
    , pPatchB(nullptr) {
    AssertLog(sdbdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundary::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundary::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundary::setPatches(Patch* patcha, Patch* patchb) {
    AssertLog(!pSetPatches);
    AssertLog(patcha != nullptr);
    AssertLog(patchb != nullptr);
    AssertLog(patcha != patchb);

    pPatchA = patcha;
    pPatchB = patchb;
    pSetPatches = true;
}

////////////////////////////////////////////////////////////////////////////////

Patch* SDiffBoundary::patchA() {
    AssertLog(pSetPatches);
    return pPatchA;
}

////////////////////////////////////////////////////////////////////////////////

Patch* SDiffBoundary::patchB() {
    AssertLog(pSetPatches);
    return pPatchB;
}

////////////////////////////////////////////////////////////////////////////////

void SDiffBoundary::setTriDirection(triangle_global_id tri, uint direction) {
    AssertLog(direction < 3);

    pTris.push_back(tri);
    pTriDirection.push_back(direction);
}

}  // namespace steps::mpi::tetopsplit
