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
#include <sstream>
#include <string>

// STEPS headers.
#include "model/spec.hpp"
#include "solver/fwd.hpp"
#include "solver/patchdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "solver/vessreacdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

VesSReacdef::VesSReacdef(Statedef* sd, vessreac_global_id idx, model::VesSReac* sr)
    : pStatedef(sd)
    , pIdx(idx)
    , pExtent(0)
    , pVessreac()
    , pOrder()
    , pKcst()
    , pImmobility()
    , pSetupdone(false)
    , pSurface_surface(true) {
    AssertLog(pStatedef != nullptr);
    AssertLog(sr != nullptr);

    pVessreac = sr;

    pName = sr->getID();
    pOrder = sr->getOrder();

    if (pOrder == 0) {
        std::ostringstream os;
        os << "\nModel contains zero-order vesicle surface reaction, which are not "
              "permitted. ";
        os << " Zero-order volume reaction may be used instead.";
        ArgErrLog(os.str());
    }

    pKcst = sr->getKcst();
    pImmobility = sr->getImmobilization();
    pMaxDistance = sr->getMaxDistance();

    pOlhs = sr->getOLHS();
    pSlhs = sr->getSLHS();
    pVlhs = sr->getVLHS();
    pLlhs = sr->getLLHS();
    pOrhs = sr->getORHS();
    pSrhs = sr->getSRHS();
    pVrhs = sr->getVRHS();
    pIrhs = sr->getIRHS();
    pLrhs = sr->getLRHS();

    pVdep = sr->getVDeps();

    // No ilhs now
    pOrient = VesSReacdef::OUTSIDE;

    uint nspecs = pStatedef->countSpecs();
    if (nspecs == 0) {  // Would be weird, but okay.
        return;
    }

    pSpec_S_LHS.container().resize(nspecs);
    pSpec_S_RHS.container().resize(nspecs);
    pSpec_S_UPD.container().resize(nspecs);

    pSpec_V_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_V_LHS.container().resize(nspecs);
    pSpec_V_RHS.container().resize(nspecs);
    pSpec_V_UPD.container().resize(nspecs);

    uint nlspecs = pStatedef->countLinkSpecs();

    pSpec_L_DEP.container().resize(nlspecs, DEP_NONE);
    pSpec_L_LHS.container().resize(nlspecs);
    pSpec_L_RHS.container().resize(nlspecs);
    pSpec_L_UPD.container().resize(nlspecs);

    pSpec_O_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_O_LHS.container().resize(nspecs);

    pSpec_O_RHS.container().resize(nspecs);
    pSpec_O_UPD.container().resize(nspecs);

    pSpec_I_RHS.container().resize(nspecs);
    pSpec_I_UPD.container().resize(nspecs);

    pSpec_VDEP.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void VesSReacdef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, pExtent);
}

////////////////////////////////////////////////////////////////////////////////

void VesSReacdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pKcst);
    util::restore(cp_file, pExtent);
}

////////////////////////////////////////////////////////////////////////////////

void VesSReacdef::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReacdef::reset() {
    pKcst = pVessreac->getKcst();
    pExtent = 0;
}

////////////////////////////////////////////////////////////////////////////////

void VesSReacdef::setup() {
    AssertLog(pSetupdone == false);

    for (auto const& ol: pOlhs) {
        pSurface_surface = false;
        AssertLog(outside());
        spec_global_id sidx = pStatedef->getSpecIdx(ol);
        pSpec_O_LHS[sidx] += 1;
    }

    for (auto const& sl: pSlhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(sl);
        pSpec_S_LHS[sidx] += 1;
    }

    for (auto const& vl: pVlhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(vl);
        pSpec_V_LHS[sidx] += 1;
    }

    for (auto const& ll: pLlhs) {
        linkspec_global_id sidx = pStatedef->getLinkSpecIdx(ll);
        pSpec_L_LHS[sidx] += 1;
    }

    for (auto const& sr: pSrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(sr);
        pSpec_S_RHS[sidx] += 1;
    }

    for (auto const& orh: pOrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(orh);
        pSpec_O_RHS[sidx] += 1;
    }

    for (auto const& vrh: pVrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(vrh);
        pSpec_V_RHS[sidx] += 1;
    }

    for (auto const& lrh: pLrhs) {
        linkspec_global_id sidx = pStatedef->getLinkSpecIdx(lrh);
        pSpec_L_RHS[sidx] += 1;
    }

    for (auto const& irh: pIrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(irh);
        pSpec_I_RHS[sidx] += 1;
    }

    for (auto const& vdep: pVdep) {
        spec_global_id sidx = pStatedef->getSpecIdx(vdep);
        pSpec_VDEP[sidx] += 1;
    }

    // Now set up the update vector
    uint ngspecs = pStatedef->countSpecs();

    // Deal with surface.
    for (auto s: spec_global_id::range(ngspecs)) {
        int lhs = static_cast<int>(pSpec_S_LHS[s]);
        int rhs = static_cast<int>(pSpec_S_RHS[s]);
        int aux = pSpec_S_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_S_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_S_UPD_Coll.push_back(s);
        }
    }

    // Deal with outside.
    for (auto s: spec_global_id::range(ngspecs)) {
        int lhs = (outside() ? static_cast<int>(pSpec_O_LHS[s]) : 0);
        int rhs = static_cast<int>(pSpec_O_RHS[s]);
        int aux = pSpec_O_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_O_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_O_UPD_Coll.push_back(s);
        }
    }

    // Deal with vesicle surface.
    for (auto s: spec_global_id::range(ngspecs)) {
        int lhs = static_cast<int>(pSpec_V_LHS[s]);
        int rhs = static_cast<int>(pSpec_V_RHS[s]);
        int aux = pSpec_V_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_V_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_V_UPD_Coll.push_back(s);
        }
    }

    // Deal with inside.
    for (auto s: spec_global_id::range(ngspecs)) {
        // No ilhs species
        uint rhs = pSpec_I_RHS[s];
        uint aux = pSpec_I_UPD[s] = rhs;
        if (aux != 0) {
            pSpec_I_UPD_Coll.push_back(s);
        }
    }

    uint nglspecs = pStatedef->countLinkSpecs();

    for (auto ls: linkspec_global_id::range(nglspecs)) {
        int lhs = static_cast<int>(pSpec_L_LHS[ls]);
        int rhs = static_cast<int>(pSpec_L_RHS[ls]);
        int aux = pSpec_L_UPD[ls] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_L_DEP[ls] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_L_UPD_Coll.push_back(ls);
        }
    }

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqOutside() const {
    AssertLog(pSetupdone == true);

    // This can be checked by seeing if DEP_O or RHS_O is non-zero
    // for any species.
    uint ngspecs = pStatedef->countSpecs();
    for (auto s: spec_global_id::range(ngspecs)) {
        if (reqspec_O(s) == true) {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqInside() const {
    AssertLog(pSetupdone == true);

    // This can be checked by seeing if DEP_I or RHS_I is non-zero
    // for any species.
    uint ngspecs = pStatedef->countSpecs();
    for (auto s: spec_global_id::range(ngspecs)) {
        if (reqspec_I(s) == true) {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqSurface() const {
    AssertLog(pSetupdone == true);

    uint ngspecs = pStatedef->countSpecs();
    for (auto s: spec_global_id::range(ngspecs)) {
        if (reqspec_S(s) == true) {
            return true;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::lhs_S(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::lhs_O(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    AssertLog(outside());
    return pSpec_O_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::lhs_V(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_V_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::lhs_L(linkspec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countLinkSpecs());
    return pSpec_L_LHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::vdep(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_VDEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::dep_S(spec_global_id gidx) {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());

    if (pSpec_S_DEP.find(gidx) == pSpec_S_DEP.end()) {
        return DEP_NONE;
    } else {
        return pSpec_S_DEP[gidx];
    }
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::dep_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    AssertLog(outside());
    return pSpec_O_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::dep_V(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_V_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::dep_L(linkspec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countLinkSpecs());
    return pSpec_L_DEP[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::rhs_S(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::rhs_O(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_O_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::rhs_V(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_V_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::rhs_I(spec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_I_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::rhs_L(linkspec_global_id gidx) const {
    AssertLog(gidx < pStatedef->countLinkSpecs());
    return pSpec_L_RHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::upd_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_S_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::upd_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_O_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::upd_V(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_V_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint VesSReacdef::upd_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    return pSpec_I_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

int VesSReacdef::upd_L(linkspec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countLinkSpecs());
    return pSpec_L_UPD[gidx];
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    // if (pSpec_S_DEP[gidx] != DEP_NONE) return true;
    if (pSpec_S_DEP.find(gidx) != pSpec_S_DEP.end()) {
        return true;
    }

    if (pSpec_S_RHS[gidx] != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqspec_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (outside()) {
        if (pSpec_O_DEP[gidx] != DEP_NONE) {
            return true;
        }
    }
    if (pSpec_O_RHS[gidx] != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqspec_V(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_V_DEP[gidx] != DEP_NONE) {
        return true;
    }
    if (pSpec_V_RHS[gidx] != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqspec_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countSpecs());
    if (pSpec_I_RHS[gidx] != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VesSReacdef::reqspec_L(linkspec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    AssertLog(gidx < pStatedef->countLinkSpecs());
    if (pSpec_L_DEP[gidx] != DEP_NONE) {
        return true;
    }
    if (pSpec_L_RHS[gidx] != 0) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
