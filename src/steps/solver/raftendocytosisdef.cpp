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

#include "solver/raftendocytosisdef.hpp"

// STL headers.
#include <string>

// STEPS headers.
#include "solver/fwd.hpp"
#include "solver/patchdef.hpp"
#include "solver/statedef.hpp"
#include "solver/types.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

RaftEndocytosisdef::RaftEndocytosisdef(Statedef* sd,
                                       raftendocytosis_global_id idx,
                                       model::RaftEndocytosis* raftendo)
    : pStatedef(sd)
    , pIdx(idx)
    , pKcst()
    , pIrhs()
    , pSetupdone(false)
    , pExtent(0)
    , pVes_I_RHS() {
    AssertLog(pStatedef != nullptr);
    AssertLog(raftendo != nullptr);

    pName = raftendo->getID();

    pKcst = raftendo->getKcst();
    pDefaultKcst = pKcst;

    pIrhs = raftendo->getRHS();

    pInner = raftendo->getInner();

    pSDeps = raftendo->getSpecDeps();

    uint nspecs = pStatedef->countSpecs();
    pSpec_S_DEP.container().resize(nspecs, DEP_NONE);
    pSpec_S_LHS.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosisdef::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, pExtent);
    util::checkpoint(cp_file, pEvents);
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosisdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pKcst);
    util::restore(cp_file, pExtent);
    util::restore(cp_file, pEvents);
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosisdef::setup() {
    AssertLog(pSetupdone == false);

    pVes_I_RHS_uint = pStatedef->getVesicleIdx(pIrhs);
    pVes_I_RHS = pStatedef->vesicledef(pVes_I_RHS_uint);

    for (auto const& sl: pSDeps) {
        spec_global_id sidx = pStatedef->getSpecIdx(sl);
        pSpec_S_LHS[sidx] += 1;
    }

    // Now set up the update vector
    // Deal with surface.
    for (auto s: spec_global_id::range(pStatedef->countSpecs())) {
        int lhs = static_cast<int>(pSpec_S_LHS[s]);
        if (lhs != 0) {
            pSpec_S_DEP[s] |= DEP_STOICH;
        }
    }

    // That's it
    pSetupdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosisdef::reset() {
    pKcst = pDefaultKcst;
    pExtent = 0;
    pEvents.clear();
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosisdef::setKcst(double kcst) {
    pKcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

Vesicledef* RaftEndocytosisdef::rhs_I_ves() const {
    AssertLog(pSetupdone == true);
    return pVes_I_RHS;
}

////////////////////////////////////////////////////////////////////////////////

vesicle_global_id RaftEndocytosisdef::rhs_I_ves_uint() const {
    AssertLog(pSetupdone == true);
    return pVes_I_RHS_uint;
}

////////////////////////////////////////////////////////////////////////////////

int RaftEndocytosisdef::dep_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    return pSpec_S_DEP.at(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool RaftEndocytosisdef::reqspec_S(spec_global_id gidx) const {
    AssertLog(pSetupdone == true);
    if (pSpec_S_DEP.at(gidx) != DEP_NONE) {
        return true;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<RaftEndocytosisEvent> RaftEndocytosisdef::getEvents() {
    std::vector<RaftEndocytosisEvent> copy(pEvents);
    pEvents.clear();
    return copy;
}

}  // namespace steps::solver
