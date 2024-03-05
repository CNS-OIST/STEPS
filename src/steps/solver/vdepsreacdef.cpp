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

#include "vdepsreacdef.hpp"

#include "model/vdepsreac.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

VDepSReacdef::VDepSReacdef(Statedef& sd, vdepsreac_global_id idx, model::VDepSReac& vdsr)
    : pStatedef(sd)
    , pIdx(idx)
    , pName(vdsr.getID())
    , pOrder(vdsr.getOrder())
    , pIlhs(vdsr.getILHS())
    , pOlhs(vdsr.getOLHS())
    , pSlhs(vdsr.getSLHS())
    , pIrhs(vdsr.getIRHS())
    , pOrhs(vdsr.getORHS())
    , pSrhs(vdsr.getSRHS())
    , pOrient(vdsr.getInner() ? VDepSReacdef::INSIDE : VDepSReacdef::OUTSIDE)
    , pVMin(vdsr._getVMin())
    , pVMax(vdsr._getVMax())
    , pDV(vdsr._getDV())
    , pVKTab(vdsr._getK()) {
    if (vdsr.getOrder() == 0) {
        std::ostringstream os;
        os << "Model contains zero-order voltage-dependent surface reaction, which "
              "are not permitted. ";
        ArgErrLog(os.str());
    }

    uint tablesize = vdsr._getTablesize();
    AssertLog(tablesize ==
              static_cast<uint>(std::floor((vdsr._getVMax() - vdsr._getVMin()) / vdsr._getDV())) +
                  1);

    uint nspecs = pStatedef.countSpecs();
    pSpec_S_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_S_LHS.container().resize(nspecs);
    if (pOrient == VDepSReacdef::INSIDE) {
        pSpec_I_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_I_LHS.container().resize(nspecs, 0);
    } else {
        pSpec_O_DEP.container().resize(nspecs, DEP_NONE);
        pSpec_O_LHS.container().resize(nspecs, 0);
    }
    pSpec_I_RHS.container().resize(nspecs);
    pSpec_S_RHS.container().resize(nspecs);
    pSpec_O_RHS.container().resize(nspecs);
    pSpec_I_UPD.container().resize(nspecs);
    pSpec_S_UPD.container().resize(nspecs);
    pSpec_O_UPD.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReacdef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReacdef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReacdef::setup() {
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
        spec_global_id sidx = pStatedef.getSpecIdx(*ol);
        pSpec_O_LHS[sidx] += 1;
    }

    for (auto const& il: pIlhs) {
        pSurface_surface = false;
        spec_global_id sidx = pStatedef.getSpecIdx(*il);
        pSpec_I_LHS[sidx] += 1;
    }

    for (auto const& sl: pSlhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*sl);
        pSpec_S_LHS[sidx] += 1;
    }

    for (auto const& ir: pIrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*ir);
        pSpec_I_RHS[sidx] += 1;
    }

    for (auto const& sr: pSrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*sr);
        pSpec_S_RHS[sidx] += 1;
    }

    for (auto const& orh: pOrhs) {
        spec_global_id sidx = pStatedef.getSpecIdx(*orh);
        pSpec_O_RHS[sidx] += 1;
    }

    // Now set up the update vector
    uint ngspecs = pStatedef.countSpecs();
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

    // This has to come before the final loop
    pSetupdone = true;

    for (auto s: spec_global_id::range(ngspecs)) {
        if (reqspec_I(s) == true) {
            pReqInside = true;
        }
        if (reqspec_O(s) == true) {
            pReqOutside = true;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double VDepSReacdef::getVDepK(double v) const {
    AssertLog(pSetupdone == true);
    AssertLog(!pVKTab.empty());
    if (v > pVMax) {
        std::ostringstream os;
        os << "Voltage to VDepSReac::getVDepRate higher than maximum: ";
        os << v << " > " << pVMax;
        ProgErrLog(os.str());
    }
    if (v < pVMin) {
        std::ostringstream os;
        os << "Voltage to VDepSReac::getVDepRate lower than minimum: ";
        os << v << " < " << pVMin;
        ProgErrLog(os.str());
    }

    double v2 = ((v - pVMin) / pDV);
    double lv = floor(v2);
    auto lvidx = static_cast<uint>(lv);
    uint uvidx = static_cast<uint>(ceil(v2));
    double r = v2 - lv;

    const auto& K = pVKTab;
    return ((1.0 - r) * K[lvidx]) + (r * K[uvidx]);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::lhs_I(spec_global_id gidx) const {
    if (outside()) {
        return 0;
    }
    return pSpec_I_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::lhs_S(spec_global_id gidx) const {
    return pSpec_S_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::lhs_O(spec_global_id gidx) const {
    if (inside()) {
        return 0;
    }
    return pSpec_O_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::dep_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (outside()) {
        return DEP_NONE;
    }
    return pSpec_I_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::dep_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::dep_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (inside()) {
        return DEP_NONE;
    }
    return pSpec_O_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::rhs_I(spec_global_id gidx) const {
    return pSpec_I_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::rhs_S(spec_global_id gidx) const {
    return pSpec_S_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint VDepSReacdef::rhs_O(spec_global_id gidx) const {
    return pSpec_O_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::upd_I(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_I_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::upd_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int VDepSReacdef::upd_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_O_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool VDepSReacdef::reqspec_I(spec_global_id gidx) const {
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

bool VDepSReacdef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_S_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return pSpec_S_RHS.at(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool VDepSReacdef::reqspec_O(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (outside()) {
        if (pSpec_O_DEP.at(gidx) != DEP_NONE) {
            return true;
        }
    }
    return pSpec_O_RHS.at(gidx) != 0;
}

}  // namespace steps::solver
