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

#include "ohmiccurrdef.hpp"

#include "model/chanstate.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

OhmicCurrdef::OhmicCurrdef(Statedef& sd, ohmiccurr_global_id gidx, model::OhmicCurr& oc)
    : pIdx(gidx)
    , pName(oc.getID())
    , pChanState(oc.getChanState().getID())
    , pG(oc.getG())
    , pERev(oc.getERev())
    , pSpec_CHANSTATE{} {
    uint nspecs = sd.countSpecs();
    pSpec_DEP.container().resize(nspecs, DEP_NONE);
    AssertLog(getG() >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurrdef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurrdef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void OhmicCurrdef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);

    spec_global_id chidx = sd.getSpecIdx(pChanState);

    pSpec_CHANSTATE = chidx;
    pSpec_DEP[chidx] |= DEP_STOICH;

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id OhmicCurrdef::chanstate() const {
    AssertLog(pSetupdone == true);
    return pSpec_CHANSTATE;
}

////////////////////////////////////////////////////////////////////////////////

int OhmicCurrdef::dep(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool OhmicCurrdef::req(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
