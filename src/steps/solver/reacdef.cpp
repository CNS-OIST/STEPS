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

#include "reacdef.hpp"

#include "compdef.hpp"
#include "model/reac.hpp"
#include "statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

Reacdef::Reacdef(Statedef& sd, reac_global_id idx, model::Reac& r)
    : pIdx(idx)
    , pName(r.getID())
    , pOrder(r.getOrder())
    , pKcst(r.getKcst())
    , pLhs(r.getLHS())
    , pRhs(r.getRHS()) {
    uint nspecs = sd.countSpecs();
    pSpec_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_LHS.container().resize(nspecs);
    pSpec_RHS.container().resize(nspecs);
    pSpec_UPD.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void Reacdef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Reacdef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void Reacdef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);

    // first copy the information about the reaction stoichiometry from Reac
    // object
    for (auto const& l: pLhs) {
        spec_global_id sidx = sd.getSpecIdx(*l);
        pSpec_LHS[sidx] += 1;
    }
    for (auto const& r: pRhs) {
        spec_global_id sidx = sd.getSpecIdx(*r);
        pSpec_RHS[sidx] += 1;
    }

    // Now set up the update vector
    for (auto i: pSpec_LHS.range()) {
        auto lhs = static_cast<int>(pSpec_LHS[i]);
        auto rhs = static_cast<int>(pSpec_RHS[i]);
        int aux = pSpec_UPD[i] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_DEP[i] |= DEP_STOICH;
        }
        if (aux != 0) {
            pSpec_UPD_Coll.emplace_back(i);
        }
    }

    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

uint Reacdef::lhs(spec_global_id gidx) const {
    return pSpec_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

depT Reacdef::dep(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint Reacdef::rhs(spec_global_id gidx) const {
    return pSpec_RHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

int Reacdef::upd(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_UPD.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool Reacdef::reqspec(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    if (pSpec_RHS.at(gidx) != 0) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
