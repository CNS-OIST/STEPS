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
#include "chandef.hpp"
#include "types.hpp"

#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

Chandef::Chandef(Statedef* sd, chan_global_id idx, model::Chan* c)
    : pStatedef(sd)
    , pIdx(idx)
    , pSetupdone(false)
    , pNChanStates(0) {
    AssertLog(pStatedef != nullptr);
    AssertLog(c != nullptr);
    pName = c->getID();

    pChanStatesVec = c->getAllChanStates();
    pNChanStates = pChanStatesVec.size();
    pChanStates.resize(pNChanStates);
}

////////////////////////////////////////////////////////////////////////////////

Chandef::~Chandef() {}

////////////////////////////////////////////////////////////////////////////////

void Chandef::checkpoint(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Chandef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Chandef::setup() {
    AssertLog(pSetupdone == false);
    AssertLog(pChanStatesVec.size() == nchanstates());
    for (uint i = 0; i < nchanstates(); ++i) {
        spec_global_id gidx = pStatedef->getSpecIdx(pChanStatesVec[i]);
        pChanStates[i] = gidx;
    }

    pSetupdone = true;
}

}  // namespace steps::solver
