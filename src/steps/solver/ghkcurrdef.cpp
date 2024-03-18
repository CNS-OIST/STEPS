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

#include "ghkcurrdef.hpp"

#include "model/chanstate.hpp"
#include "model/ghkcurr.hpp"
#include "model/spec.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

GHKcurrdef::GHKcurrdef(Statedef& sd, ghkcurr_global_id gidx, model::GHKcurr& ghk)
    : pIdx(gidx)
    , pName(ghk.getID())
    , pChanState(ghk.getChanState().getID())
    , pIon(ghk.getIon().getID())
    , pRealFlux(ghk._realflux())
    , pVirtual_oconc(ghk._voconc())
    , pVshift(ghk._vshift())
    , pPerm(ghk._P())
    , pValence(ghk._valence()) {
    if (!ghk._infosupplied()) {
        std::ostringstream os;
        os << "\nPermeability not defined for GHK current object.";
        ArgErrLog(os.str());
    }

    AssertLog(pValence != 0);
    AssertLog(pPerm > 0.0);

    uint nspecs = sd.countSpecs();
    pSpec_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_VOL_DEP.container().resize(nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::checkpoint(std::fstream& /*cp_file*/) const {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurrdef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);

    spec_global_id chidx = sd.getSpecIdx(pChanState);
    spec_global_id ionidx = sd.getSpecIdx(pIon);

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
