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

#include "mpi/tetvesicle/diffboundary.hpp"

// STEPS headers.
#include "mpi/tetvesicle/comp_rdef.hpp"
#include "solver/diffboundarydef.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

DiffBoundary::DiffBoundary(solver::DiffBoundarydef* dbdef)
    : pDiffBoundarydef(dbdef)
    , pSetComps(false)
    , pCompA(nullptr)
    , pCompB(nullptr) {
    AssertLog(dbdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

DiffBoundary::~DiffBoundary() = default;

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::setComps(CompRDEF* compa, CompRDEF* compb) {
    AssertLog(pSetComps == false);
    AssertLog(compa != nullptr);
    AssertLog(compb != nullptr);
    AssertLog(compa != compb);

    pCompA = compa;
    pCompB = compb;
    pSetComps = true;
}

////////////////////////////////////////////////////////////////////////////////

CompRDEF* DiffBoundary::compA() {
    AssertLog(pSetComps == true);
    return pCompA;
}

////////////////////////////////////////////////////////////////////////////////

CompRDEF* DiffBoundary::compB() {
    AssertLog(pSetComps == true);
    return pCompB;
}

////////////////////////////////////////////////////////////////////////////////

void DiffBoundary::setTetDirection(tetrahedron_global_id tet, uint direction) {
    AssertLog(direction < 4);

    pTets.push_back(tet);
    pTetDirection.push_back(direction);
}

}  // namespace steps::mpi::tetvesicle
