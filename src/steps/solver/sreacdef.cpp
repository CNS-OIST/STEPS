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
#include "sreacdef.hpp"

#include "geom/patch.hpp"
#include "model/sreac.hpp"
#include "patchdef.hpp"
#include "statedef.hpp"
#include "types.hpp"
#include "util/error.hpp"
// logging
#include <easylogging++.h>

namespace steps::solver {

SReacdef::SReacdef(Statedef* sd, sreac_global_id idx, model::SReac* sr)
    : pStatedef(sd)
    , pIdx(idx)
    , pOrder()
    , pKcst()
    , pSetupdone(false)
    , pSurface_surface(true) {
    AssertLog(pStatedef != nullptr);
    AssertLog(sr != nullptr);

    pName = sr->getID();
    pOrder = sr->getOrder();

    if (pOrder == 0) {
        std::ostringstream os;
        os << "\nModel contains zero-order surface reaction, which are not "
              "permitted. ";
        os << " Zero-order volume reaction may be used instead\n.";
        ArgErrLog(os.str());
    }

    pKcst = sr->getKcst();
    pIlhs = sr->getILHS();
    pOlhs = sr->getOLHS();
    pSlhs = sr->getSLHS();

    pIrhs = sr->getIRHS();
    pOrhs = sr->getORHS();
    pSrhs = sr->getSRHS();

    if (sr->getInner() == true) {
        pOrient = SReacdef::INSIDE;
    } else {
        pOrient = SReacdef::OUTSIDE;
    }

    uint nspecs = pStatedef->countSpecs();

    pSpec_S_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_S_LHS.container().resize(nspecs);
    if (pOrient == SReacdef::INSIDE) {
        pSpec_I_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_I_LHS.container().resize(nspecs);
    } else {
        pSpec_O_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_O_LHS.container().resize(nspecs);
    }
    pSpec_I_RHS.container().resize(nspecs);
    pSpec_S_RHS.container().resize(nspecs);
    pSpec_O_RHS.container().resize(nspecs);
    pSpec_I_UPD.container().resize(nspecs);
    pSpec_S_UPD.container().resize(nspecs);
    pSpec_O_UPD.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void SReacdef::checkpoint(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

void SReacdef::restore(std::fstream& /*cp_file*/) {
    // reserve
}

////////////////////////////////////////////////////////////////////////////////

std::string const SReacdef::name() const {
    return pName;
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::order() const {
    return pOrder;
}

////////////////////////////////////////////////////////////////////////////////

double SReacdef::kcst() const {
    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

void SReacdef::setup() {
    AssertLog(pSetupdone == false);

    if (outside()) {
        AssertLog(pIlhs.empty());
    } else if (inside()) {
        AssertLog(pOlhs.empty());
    } else {
        AssertLog(false);
    }

    for (auto const& ol: pOlhs) {
        pSurface_surface = false;
        spec_global_id sidx = pStatedef->getSpecIdx(ol);
        pSpec_O_LHS[sidx] += 1;
    }

    for (auto const& il: pIlhs) {
        pSurface_surface = false;
        spec_global_id sidx = pStatedef->getSpecIdx(il);
        pSpec_I_LHS[sidx] += 1;
    }

    for (auto const& sl: pSlhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(sl);
        pSpec_S_LHS[sidx] += 1;
    }

    for (auto const& ir: pIrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(ir);
        pSpec_I_RHS[sidx] += 1;
    }

    for (auto const& sr: pSrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(sr);
        pSpec_S_RHS[sidx] += 1;
    }

    for (auto const& orh: pOrhs) {
        spec_global_id sidx = pStatedef->getSpecIdx(orh);
        pSpec_O_RHS[sidx] += 1;
    }

    // Now set up the update vector
    uint ngspecs = pStatedef->countSpecs();

    // Deal with surface.
    for (auto s: spec_global_id::range(ngspecs)) {
        auto lhs = static_cast<int>(pSpec_S_LHS[s]);
        auto rhs = static_cast<int>(pSpec_S_RHS[s]);
        int aux = pSpec_S_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_S_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_S_UPD_Coll.push_back(s);
        }
    }

    // Deal with inside.
    for (auto s: spec_global_id::range(ngspecs)) {
        int lhs = (inside() ? static_cast<int>(pSpec_I_LHS[s]) : 0);
        auto rhs = static_cast<int>(pSpec_I_RHS[s]);
        int aux = pSpec_I_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_I_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_I_UPD_Coll.push_back(s);
        }
    }

    // Deal with outside.
    for (auto s: spec_global_id::range(ngspecs)) {
        int lhs = (outside() ? static_cast<int>(pSpec_O_LHS[s]) : 0);
        auto rhs = static_cast<int>(pSpec_O_RHS[s]);
        int aux = pSpec_O_UPD[s] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_O_DEP[s] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_O_UPD_Coll.push_back(s);
        }
    }

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacdef::reqInside() const {
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

bool SReacdef::reqOutside() const {
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

uint SReacdef::lhs_I(spec_global_id gidx) const {
    if (outside()) {
        return 0;
    }
    return pSpec_I_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::lhs_S(spec_global_id gidx) const {
    return pSpec_S_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::lhs_O(spec_global_id gidx) const {
    if (inside()) {
        return 0;
    }
    return pSpec_O_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::dep_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (outside()) {
        return DEP_NONE;
    }
    return pSpec_I_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::dep_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::dep_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (inside()) {
        return DEP_NONE;
    }
    return pSpec_O_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::rhs_I(spec_global_id gidx) const {
    return pSpec_I_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::rhs_S(spec_global_id gidx) const {
    return pSpec_S_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint SReacdef::rhs_O(spec_global_id gidx) const {
    return pSpec_O_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::upd_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_I_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::upd_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int SReacdef::upd_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_O_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool SReacdef::reqspec_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (inside()) {
        if (pSpec_I_DEP.at(gidx) != DEP_NONE) {
            return true;
        }
    }
    if (pSpec_I_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacdef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_S_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    if (pSpec_S_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacdef::reqspec_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (outside()) {
        if (pSpec_O_DEP.at(gidx) != DEP_NONE) {
            return true;
        }
    }
    if (pSpec_O_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
