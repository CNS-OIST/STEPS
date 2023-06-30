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
#include <cassert>
#include <string>

// STEPS headers.
#include "model/linkspec.hpp"
#include "solver/linkspecdef.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
// #include "solver/statedef.hpp"
#include "solver/types.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

LinkSpecdef::LinkSpecdef(Statedef* sd, linkspec_global_id idx, model::LinkSpec* l)
    : pStatedef(sd)
    , pIdx(idx)
    , pSetupdone(false) {
    AssertLog(pStatedef != nullptr);
    AssertLog(l != nullptr);
    pName = l->getID();

    pDcst = l->getDcst();
}


////////////////////////////////////////////////////////////////////////////////

LinkSpecdef::~LinkSpecdef() = default;

////////////////////////////////////////////////////////////////////////////////

void LinkSpecdef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpecdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void LinkSpecdef::setup() {
    pSetupdone = true;
}

}  // namespace steps::solver
