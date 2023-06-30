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

#include "mpi/tetvesicle/raftendocytosis.hpp"

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/comp_vesraft.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

RaftEndocytosis::RaftEndocytosis(solver::RaftEndocytosisdef* endodef, Raft* raft)
    : pRaftEndocytosisdef(endodef)
    , pRaft(raft)
    , pCcst(0.0)
    , pKcst(0.0)
    , pActive(true)
    , pFails(0) {
    AssertLog(pRaftEndocytosisdef != nullptr);
    AssertLog(pRaft != nullptr);

    double kcst = pRaftEndocytosisdef->kcst();
    AssertLog(kcst >= 0.0);
    pKcst = kcst;

    // It's a type of 1st order reaction
    pCcst = kcst;
}

RaftEndocytosis::RaftEndocytosis(solver::RaftEndocytosisdef* endodef,
                                 Raft* raft,
                                 std::fstream& cp_file)
    : pRaftEndocytosisdef(endodef)
    , pRaft(raft) {
    AssertLog(pRaftEndocytosisdef != nullptr);
    AssertLog(pRaft != nullptr);

    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    util::restore(cp_file, pActive);
    util::restore(cp_file, pFails);
}

////////////////////////////////////////////////////////////////////////////////

RaftEndocytosis::~RaftEndocytosis() = default;

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    util::checkpoint(cp_file, pActive);
    util::checkpoint(cp_file, pFails);
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::reset() {
    resetCcst();
    setActive(true);

    pFails = 0;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::resetCcst() {
    double kcst = pRaftEndocytosisdef->kcst();
    AssertLog(kcst >= 0.0);
    pKcst = kcst;

    // It's always a 1st order reaction
    pCcst = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
    // It's always a 1st order reaction
    pCcst = k;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::setActive(bool active) {
    pActive = active;
}

////////////////////////////////////////////////////////////////////////////////

double RaftEndocytosis::rate() {
    if (inactive()) {
        return 0.0;
    }
    // I know that this uses global indices

    const auto& lhs_s_vec = endodef()->lhs_S();

    const auto& cnt_vec = pRaft->pools_global();

    for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
        uint lhs = lhs_s_vec[sg];
        if (lhs == 0) {
            continue;
        }
        uint cnt = cnt_vec[sg];

        // We have summed all available species from each tri. Compare to required
        // lhs
        if (lhs > cnt) {
            //  The required species are not available
            return 0.0;
        }
    }

    // If we got here all species are available
    // Special kind of reaction that is like a pseudo-first order reaction
    return pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void RaftEndocytosis::apply() {
    // First see if we can do the raftendocytosis. If not just exit.
    solver::Vesicledef* ves_def = endodef()->rhs_I_ves();
    solver::vesicle_global_id ves_idx = endodef()->rhs_I_ves_uint();

    TriVesRaft* tri = pRaft->tri_central();

    // Get the patch and comp pointers but check consistency
    TetVesRaft* tet;
    CompVesRaft* comp;
    if (inner()) {
        tet = tri->iTet();
    } else {
        tet = tri->oTet();
    }
    comp = tet->getCompVesRaft();

    // Store map of species lidx in patch to population in triangles
    std::map<solver::spec_global_id, int> specs_vesicle;

    // Find the position
    math::point3d pos;

    const math::point3d& tet_baryc = tet->position();
    const math::point3d& tri_baryc = tri->position();
    math::point3d tri_norm = tri->norm();

    solver::vesicle_global_id ves_gidx = pRaftEndocytosisdef->rhs_I_ves_uint();
    solver::Vesicledef* vesdef = pRaftEndocytosisdef->statedef()->vesicledef(ves_gidx);
    double ves_diameter = vesdef->diameter();
    math::point3d baryc_baryc_vec = tet_baryc - tri_baryc;
    double dotproduct = math::dot(baryc_baryc_vec, tri_norm);
    if (dotproduct < 0.0) {
        tri_norm[0] = -tri_norm[0];
        tri_norm[1] = -tri_norm[1];
        tri_norm[2] = -tri_norm[2];
    }

    // Tri normal is already normalised as returned (mg 1)
    pos = (tri_baryc + tri_norm * (ves_diameter / 2.0));

    // First we need to try the endocytosis:
    solver::vesicle_individual_id ves_unique_index;

    math::position_abs pos_abs{pos};

    ves_unique_index = comp->addVesicle(ves_def, pos_abs);

    if (ves_unique_index.valid()) {
        // I know that this uses global indices
        uint nspecs_g = endodef()->countSpecs_gobal();

        for (auto spec_gidx: solver::spec_global_id::range(nspecs_g)) {
            uint nmolcs = pRaft->pools_global()[spec_gidx];

            if (nmolcs > 0) {
                // vesicles store global species indices
                specs_vesicle[spec_gidx] += nmolcs;
            }
            // don't need to do anything else because raft will be destroyed
        }

        comp->addVesicleSpecs(ves_idx, ves_unique_index, specs_vesicle);

        // Now this is controlled directly by Raft, all updates and
        // everything are controlled there, but raft needs to know if this
        // went ahead or not.
        pRaft->appliedEndo();

        endodef()->addEvent(endodef()->statedef()->time(),
                            pRaft->getUniqueIndex(),
                            tri->idx(),
                            ves_unique_index);
    }

    else {
        pFails += 1;
        CLOG(WARNING, "general_log")
            << "\nRaftEndocytosis " << endodef()->name() << " failed. " << pFails << " Fails.\n";
    }
}

}  // namespace steps::mpi::tetvesicle
