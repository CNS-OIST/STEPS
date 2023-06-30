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
#include "ghkcurrdef.hpp"
#include "math/ghk.hpp"
#include "model/chanstate.hpp"
#include "model/spec.hpp"
#include "util/error.hpp"

// logging
#include <easylogging++.h>

namespace steps::solver {

GHKcurrdef::GHKcurrdef(Statedef* sd, ghkcurr_global_id gidx, model::GHKcurr* ghk)
    : pStatedef(sd)
    , pIdx(gidx) {
    AssertLog(pStatedef != nullptr);
    AssertLog(ghk != nullptr);

    pName = ghk->getID();
    pChanState = ghk->getChanState()->getID();
    pIon = ghk->getIon()->getID();
    pRealFlux = ghk->_realflux();
    pVirtual_oconc = ghk->_voconc();
    pVshift = ghk->_vshift();

    if (!ghk->_infosupplied()) {
        std::ostringstream os;
        os << "\nPermeability not defined for GHK current object.";
        ArgErrLog(os.str());
    }

    pValence = ghk->_valence();
    AssertLog(pValence != 0);

    double perm = ghk->_P();
    AssertLog(perm > 0.0);
    pPerm = perm;

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) {
        return;  // Would be weird, but okay.
    }

    pSpec_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_VOL_DEP.container().resize(nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::setup() {
    AssertLog(pSetupdone == false);

    spec_global_id chidx = pStatedef->getSpecIdx(pChanState);
    spec_global_id ionidx = pStatedef->getSpecIdx(pIon);

    pSpec_CHANSTATE = chidx;
    pSpec_ION = ionidx;

    pSpec_DEP[chidx] |= DEP_STOICH;
    // Note: The dependency here is indirect. The flux is not modelled as
    // a 2nd order reaction ion + channel -> movement of ion
    // instead the flux comes from the GHK flux equation and is then
    // modelled as a 1st order reaction, so the channel is the only dependency
    // Update is a difference matter though, each event will change the number of
    // molecules in bordering tetrahedrons.
    // And, the concentration of ion in the inner and outer compartment
    // has an affect on the RATE

    // And, the concentration of ion in the inner and outer compartment
    // has an affect on the RATE
    pSpec_VOL_DEP[ionidx] |= DEP_RATE;

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id GHKcurrdef::chanstate() const {
    AssertLog(pSetupdone == true);
    return pSpec_CHANSTATE;
}

////////////////////////////////////////////////////////////////////////////////

spec_global_id GHKcurrdef::ion() const {
    AssertLog(pSetupdone == true);
    return pSpec_ION;
}
////////////////////////////////////////////////////////////////////////////////

int GHKcurrdef::dep(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int GHKcurrdef::dep_v(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_VOL_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool GHKcurrdef::req(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool GHKcurrdef::req_v(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_VOL_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
