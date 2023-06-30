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
#include "solver/exocytosisdef.hpp"
#include "solver/fwd.hpp"
#include "solver/patchdef.hpp"
#include "solver/statedef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

Exocytosisdef::Exocytosisdef(Statedef* sd, exocytosis_global_id idx, model::Exocytosis* exo)
    : pStatedef(sd)
    , pIdx(idx)
    , pKcst()
    , pExocytosis(exo)
    , pSetupdone(false)
    , pExtent(0)
    , pRaft()
    , pRaftdef(nullptr) {
    AssertLog(pStatedef != nullptr);
    AssertLog(exo != nullptr);

    pName = exo->getID();

    pKcst = exo->getKcst();

    pVDeps = exo->getSpecDeps();

    pRaft = exo->getRaft();

    pKissAndRun = exo->getKissAndRun();

    for (auto const& specs: exo->getKissAndRunSpecChanges()) {
        spec_global_id sidx_src = pStatedef->getSpecIdx(specs.first);
        spec_global_id sidx_dst = pStatedef->getSpecIdx(specs.second);
        AssertLog(sidx_src != sidx_dst);
        pKissAndRunSpecChanges[sidx_src] = sidx_dst;
    }

    if (pRaft != nullptr) {
        AssertLog(pKissAndRun == false);
    }

    uint nspecs = statedef()->countSpecs();
    if (nspecs == 0) {
        return;  // Would be weird, but okay.
    }

    pSpec_V_LHS.container().resize(nspecs);

    pSpec_V_DEP.container().resize(nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

Exocytosisdef::~Exocytosisdef() {}

////////////////////////////////////////////////////////////////////////////////

void Exocytosisdef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, pExtent);
    util::checkpoint(cp_file, pEvents);
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosisdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pKcst);
    util::restore(cp_file, pExtent);
    util::restore(cp_file, pEvents);
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosisdef::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosisdef::setup() {
    AssertLog(pSetupdone == false);

    for (auto const& sl: pVDeps) {
        spec_global_id sidx = pStatedef->getSpecIdx(sl);
        pSpec_V_LHS[sidx] += 1;
    }

    // Now set up the update vector
    uint nspecs = pStatedef->countSpecs();
    // Deal with surface.
    // most convenient just to work with uints here (not strong ids)
    for (auto i: spec_global_id::range(nspecs)) {
        int lhs = static_cast<int>(pSpec_V_LHS[i]);
        // int rhs = static_cast<int>(pSpec_S_RHS[i]);
        // int aux = pSpec_S_UPD[i] = (rhs - lhs);
        if (lhs != 0) {
            pSpec_V_DEP[i] |= DEP_STOICH;
        }
        // if (aux != 0) pSpec_S_UPD_Coll.push_back(i);
    }

    // Find the raftdef. Might be NULL
    if (pRaft != nullptr) {
        raft_global_id raftidx = pStatedef->getRaftIdx(pRaft);
        pRaftdef = pStatedef->rafts()[raftidx];
    }

    // That's it
    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

int Exocytosisdef::dep_V(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_V_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

uint Exocytosisdef::lhs_V(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_V_LHS.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<ExocytosisEvent> Exocytosisdef::getEvents() {
    std::vector<ExocytosisEvent> copy(pEvents);
    pEvents.clear();
    return copy;
}

}  // namespace steps::solver
