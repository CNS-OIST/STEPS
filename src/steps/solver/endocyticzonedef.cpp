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

// STL headers.
#include <iostream>
#include <string>

// STEPS headers.
#include "solver/endocyticzonedef.hpp"
#include "solver/fwd.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

namespace steps::solver {

EndocyticZonedef::EndocyticZonedef(Statedef* sd, tetmesh::EndocyticZone* z)
    : pStatedef(sd) {
    AssertLog(pStatedef != nullptr);
    AssertLog(z != nullptr);

    pName = z->getID();

    pTris = z->getAllTriIndices();
}

////////////////////////////////////////////////////////////////////////////////

void EndocyticZonedef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void EndocyticZonedef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

EndocyticZonedef::~EndocyticZonedef() = default;

}  // namespace steps::solver
