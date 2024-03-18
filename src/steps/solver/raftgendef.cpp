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

#include "solver/raftgendef.hpp"

#include "solver/statedef.hpp"
#include "util/error.hpp"

namespace steps::solver {

RaftGendef::RaftGendef(Statedef& sd, raftgen_global_id idx, model::RaftGen& raftgen)
    : pIdx(idx)
    , pName(raftgen.getID())
    , pKcst(raftgen.getKcst())
    , pCountSpecs(sd.countSpecs())
    , pSDeps(raftgen.getSpecSignature())
    , pRaft(raftgen.getRaft()) {
    pSpec_S_DEP.container().resize(pCountSpecs, DEP_NONE);
    pSpec_S_LHS.container().resize(pCountSpecs);
}


////////////////////////////////////////////////////////////////////////////////

void RaftGendef::checkpoint(std::fstream& /*cp_file*/) const {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void RaftGendef::restore(std::fstream& /*cp_file*/) {
    // Reserve
}

////////////////////////////////////////////////////////////////////////////////

void RaftGendef::setup(const Statedef& sd) {
    AssertLog(pSetupdone == false);

    for (auto const& sl: pSDeps) {
        spec_global_id sidx = sd.getSpecIdx(*sl);
        pSpec_S_LHS[sidx] += 1;
    }

    // Now set up the update vector
    // Deal with surface.
    for (auto s: spec_global_id::range(sd.countSpecs())) {
        int lhs = static_cast<int>(pSpec_S_LHS[s]);
        if (lhs != 0) {
            pSpec_S_DEP[s] |= DEP_STOICH;
        }
    }

    // Find the raftdef.
    raft_global_id raftidx = sd.getRaftIdx(pRaft);
    pRaftdef = sd.rafts()[raftidx].get();

    // That's it
    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

int RaftGendef::dep_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool RaftGendef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_S_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return false;
}

}  // namespace steps::solver
