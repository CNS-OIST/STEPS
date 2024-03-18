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

#include "chandef.hpp"

#include "model/chan.hpp"
#include "model/chanstate.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

Chandef::Chandef(Statedef& /*sd*/, chan_global_id idx, model::Chan& c)
    : pIdx(idx)
    , pName(c.getID())
    , pChanStatesVec(c.getAllChanStates()) {
    pChanStates.resize(pChanStatesVec.size());
}

////////////////////////////////////////////////////////////////////////////////

void Chandef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Chandef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Chandef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);
    const auto& chan_states = pChanStatesVec;
    AssertLog(chan_states.size() == nchanstates());
    for (uint i = 0; i < nchanstates(); ++i) {
        spec_global_id gidx = sd.getSpecIdx(*chan_states[i]);
        pChanStates[i] = gidx;
    }

    pSetupdone = true;
}

}  // namespace steps::solver
