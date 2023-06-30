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
#include "tet.hpp"
#include "tetode.hpp"
#include "tri.hpp"
#include "util/checkpointing.hpp"

namespace steps::tetode {

Tet::Tet(tetrahedron_global_id idx,
         solver::Compdef* cdef,
         double vol,
         double a0,
         double a1,
         double a2,
         double a3,
         double d0,
         double d1,
         double d2,
         double d3,
         tetrahedron_global_id tet0,
         tetrahedron_global_id tet1,
         tetrahedron_global_id tet2,
         tetrahedron_global_id tet3)
    : pCompdef(cdef)
    , pIdx(idx)
    , pVol(vol)
    , pTets()
    , pNextTri()
    , pNextTet()
    , pAreas()
    , pDist() {
    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);


    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i = 0; i <= 3; ++i) {
        pNextTet[i] = nullptr;
        pNextTri[i] = nullptr;
    }
    pTets[0] = tet0;
    pTets[1] = tet1;
    pTets[2] = tet2;
    pTets[3] = tet3;

    pAreas[0] = a0;
    pAreas[1] = a1;
    pAreas[2] = a2;
    pAreas[3] = a3;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
    pDist[3] = d3;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTet(uint i, Tet* t) {
    if (t->compdef() != compdef()) {
        pNextTet[i] = nullptr;
    } else {
        pNextTet[i] = t;
        if (pNextTri[i] != nullptr) {
            CLOG(INFO, "general_log") << "WARNING: writing over nextTri index " << i;
        }
        pNextTri[i] = nullptr;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(uint i, Tri* t) {
    pNextTet[i] = nullptr;
    pNextTri[i] = t;
}


////////////////////////////////////////////////////////////////////////////////

void Tet::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void Tet::restore(std::fstream& /*cp_file*/) {
    // reserve
}

}  // namespace steps::tetode
