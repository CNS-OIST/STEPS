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

#include "mpi/tetvesicle/raftdis.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

RaftDis::RaftDis(solver::RaftDisdef* rgdef, Raft* raft)
    : pRaftDisdef(rgdef)
    , pRaft(raft)
    , rExtent(0)
    , pActive(true) {
    AssertLog(pRaftDisdef != nullptr);

    double kcst = pRaftDisdef->kcst();
    AssertLog(kcst >= 0.0);
    pKcst = kcst;
}

RaftDis::RaftDis(solver::RaftDisdef* rgdef, Raft* raft, std::fstream& cp_file)
    : pRaftDisdef(rgdef)
    , pRaft(raft) {
    AssertLog(pRaftDisdef != nullptr);

    util::restore(cp_file, pKcst);
    util::restore(cp_file, rExtent);
    util::restore(cp_file, pActive);
}

////////////////////////////////////////////////////////////////////////////////

RaftDis::~RaftDis() = default;

////////////////////////////////////////////////////////////////////////////////

void RaftDis::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, rExtent);
    util::checkpoint(cp_file, pActive);
}


////////////////////////////////////////////////////////////////////////////////

void RaftDis::reset() {
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void RaftDis::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
}

////////////////////////////////////////////////////////////////////////////////

double RaftDis::rate() {
    if (inactive()) {
        return 0.0;
    }
    // This uses global indices already

    const auto& lhs_s_vec = def()->lhs_S();

    for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
        uint lhs = lhs_s_vec[sg];
        if (lhs == 0) {
            continue;
        }
        uint cnt = pRaft->pools_global()[sg];

        //  Compare to required lhs
        if (cnt > lhs) {
            //  Above the threshold
            return 0.0;
        }
    }

    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

void RaftDis::apply() {
    TriVesRaft* tri = raft()->tri_central();

    // Does a global update so not setting connected tets right now

    solver::Patchdef* pdef = tri->patchdef();

    for (auto spec_gidx: solver::spec_global_id::range(def()->countSpecs_global())) {
        uint nmolcs = pRaft->pools_global()[spec_gidx];

        if (nmolcs > 0) {
            // Set the tri count
            solver::spec_local_id spec_lidx = pdef->specG2L(spec_gidx);

            if (spec_lidx.unknown()) {
                std::ostringstream os;
                os << "\nCan't apply RaftDis: species undefined in patch.";
                ArgErrLog(os.str());
            }

            uint prev_count = tri->pools()[spec_lidx];

            uint new_count = prev_count + nmolcs;

            AssertLog(new_count > 0);

            // Now done on vesicle dts so no need for period argument
            tri->setCount(spec_lidx, new_count);
        }
    }

    // nothing else to do because patch now takes care of all the cleanup
}

}  // namespace steps::mpi::tetvesicle
