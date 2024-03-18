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

#include "mpi/tetvesicle/raft.hpp"

// STEPS headers.
#include "math/point.hpp"
#include "mpi/tetvesicle/patch_vesraft.hpp"
#include "mpi/tetvesicle/raftdis.hpp"
#include "mpi/tetvesicle/raftendocytosis.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "solver/raftdef.hpp"
#include "util/checkpointing.hpp"
#include "util/vocabulary.hpp"

namespace steps::mpi::tetvesicle {

Raft::Raft(solver::Raftdef* raftdef,
           PatchVesRaft* patch,
           TriVesRaft* central_tri,
           math::position_abs& pos,
           solver::raft_individual_id unique_index)
    : pDef(raftdef)
    , pPatch(patch)
    , pIndex(unique_index)
    , pPos(pos)
    , pTri_central(central_tri)
    , pScaledDcst(0.0)
    , pAppliedEndo(false)
    , pImmobility(0) {
    AssertLog(pDef != nullptr);

    RaftDiss.container().resize(def()->countRaftDiss());
    RaftEndocytosiss.container().resize(def()->countRaftEndocytosis());
    pPoolCount.container().resize(def()->countSpecs_global());

    // Now time to setup the diffusion stuff

    TriVesRaft* next[3] = {
        pTri_central->nextTri(0),
        pTri_central->nextTri(1),
        pTri_central->nextTri(2),
    };

    std::array<double, 3> d = {0.0, 0.0, 0.0};
    for (uint i = 0; i < d.size(); ++i) {
        // Compute the scaled diffusion constant.
        double dist = pTri_central->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if (next[i]->patchdef() == pTri_central->patchdef()) {
                d[i] = (pTri_central->length(i) * getDcst()) / (pTri_central->area() * dist);
            } else {
                d[i] = 0.0;
            }
        }
    }

    // Compute scaled "diffusion constant".
    for (uint i = 0; i < d.size(); ++i) {
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    } else {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }

    // Set the overlap
    pTris_overlap_vec = pPatch->getRaftOverlap(def()->gidx(), pTri_central->idx());

    pPatch->solverVesRaft()->recordRaft_(pIndex, this);
}

////////////////////////////////////////////////////////////////////////////////

Raft::Raft(solver::Raftdef* raftdef, PatchVesRaft* patch, std::fstream& cp_file)
    : pDef(raftdef)
    , pPatch(patch)
    , pAppliedEndo(false) {
    RaftDiss.container().resize(def()->countRaftDiss());
    RaftEndocytosiss.container().resize(def()->countRaftEndocytosis());

    util::restore(cp_file, pIndex);
    util::restore(cp_file, pPos);
    triangle_global_id tri_central;
    util::restore(cp_file, tri_central);
    pTri_central = patch->solverVesRaft()->tri_(tri_central);
    util::restore(cp_file, pTris_overlap_vec);
    util::restore(cp_file, pPoolCount);
    util::restore(cp_file, pScaledDcst);
    util::restore(cp_file, pCDFSelector);
    util::restore(cp_file, pImmobility);
    util::restore(cp_file, pRaftSReac_inactive);

    pPatch->solverVesRaft()->recordRaft_(pIndex, this);

    // Create endocytosis kprocs. These belong to raft because they are associated
    // with stuff in the raft surface
    uint nendos = def()->countRaftEndocytosis();
    for (auto i: solver::raftendocytosis_local_id::range(nendos)) {
        auto& endodef = def()->raftendocytosisdef(i);
        auto* endo = new RaftEndocytosis(&endodef, this, cp_file);
        raftendos()[i] = endo;
    }

    // Create raft dis kprocs
    uint ndiss = def()->countRaftDiss();
    for (auto i: solver::raftdis_local_id::range(ndiss)) {
        auto& rddef = def()->raftdisdef(i);
        auto* rdis = new RaftDis(&rddef, this, cp_file);
        raftdiss()[i] = rdis;
    }
}

////////////////////////////////////////////////////////////////////////////////

Raft::~Raft() {
    // Now raft does not get destroyed during raftdis nor raftendo application.
    // Instead, the patch checks for application then goes ahead and
    // takes care of all the necessary cleanup
    for (auto rd: RaftDiss) {
        delete rd;
    }

    for (auto re: RaftEndocytosiss) {
        delete re;
    }

    // Need to remove overlap from compartment
    for (auto const& t: getOverlapVec()) {
        // t is global index, patch now needs local index
        triangle_local_id tlidx = pPatch->triidx_G_to_L(t);

        pPatch->tri(tlidx)->removeRaftref(this);
    }

    patch()->solverVesRaft()->removeRaft_(getUniqueIndex(), this);
}

////////////////////////////////////////////////////////////////////////////////

void Raft::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pIndex);
    util::checkpoint(cp_file, pPos);
    util::checkpoint(cp_file, pTri_central->idx());
    util::checkpoint(cp_file, pTris_overlap_vec);
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pScaledDcst);
    util::checkpoint(cp_file, pCDFSelector);
    util::checkpoint(cp_file, pImmobility);
    util::checkpoint(cp_file, pRaftSReac_inactive);
    // NOTE only time pAppliedEndo can be true is when Raft is scheduled for deletion, which should
    // take place before the call to checkpoint

    for (auto const& re: RaftEndocytosiss) {
        re->checkpoint(cp_file);
    }

    for (auto const& rd: RaftDiss) {
        rd->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Raft::setupKProcs() {
    // Create endocytosis kprocs. These belong to raft because they are associated
    // with stuff in the raft surface
    uint nendos = def()->countRaftEndocytosis();
    for (auto i: solver::raftendocytosis_local_id::range(nendos)) {
        auto& endodef = def()->raftendocytosisdef(i);
        auto* endo = new RaftEndocytosis(&endodef, this);
        raftendos()[i] = endo;
    }

    // Create raft dis kprocs
    uint ndiss = def()->countRaftDiss();
    for (auto i: solver::raftdis_local_id::range(ndiss)) {
        auto& rddef = def()->raftdisdef(i);
        auto* rdis = new RaftDis(&rddef, this);
        raftdiss()[i] = rdis;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Raft::setPosition(math::position_abs& new_pos, TriVesRaft* tri_central) {
    // TODO Error checking?
    pPos = new_pos;
    pTri_central = tri_central;

    // First reset the dcst
    pScaledDcst = 0.0;

    pTris_overlap_vec = patch()->getRaftOverlap(def()->gidx(), tri_central->idx());

    // Now time to setup the diffusion stuff

    TriVesRaft* next[3] = {
        pTri_central->nextTri(0),
        pTri_central->nextTri(1),
        pTri_central->nextTri(2),
    };

    double d[3] = {0.0, 0.0, 0.0};
    for (uint i = 0; i < 3; ++i) {
        // Compute the scaled diffusion constant.
        double dist = pTri_central->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if (next[i]->patchdef() == pTri_central->patchdef()) {
                d[i] = (pTri_central->length(i) * getDcst()) / (pTri_central->area() * dist);
            } else {
                d[i] = 0.0;
            }
        }
    }

    // Compute scaled "diffusion constant".
    for (uint i = 0; i < 3; ++i) {
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    } else {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

TriVesRaft* Raft::selectDirectionTri(double unf) {
    int iSel = 0;
    for (; iSel < 2; ++iSel) {
        if (unf < pCDFSelector[iSel]) {
            break;
        }
    }

    // Direction iSel.
    TriVesRaft* nexttri = pTri_central->nextTri(iSel);

    return nexttri;
}

////////////////////////////////////////////////////////////////////////////////

void Raft::setSpecCountByLidx(solver::spec_local_id slidx, uint count) {
    AssertLog(slidx < def()->countSpecs());
    solver::spec_global_id spec_gidx = def()->specL2G(slidx);
    pPoolCount[spec_gidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

uint Raft::getSpecCountByLidx(solver::spec_local_id slidx) {
    AssertLog(slidx < def()->countSpecs());
    solver::spec_global_id spec_gidx = def()->specL2G(slidx);
    return pPoolCount[spec_gidx];
}

////////////////////////////////////////////////////////////////////////////////

void Raft::setSpecCountByGidx(solver::spec_global_id sgidx, uint count) {
    AssertLog(sgidx < def()->countSpecs_global());
    pPoolCount[sgidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void Raft::incSpecCountByGidx(solver::spec_global_id sgidx, uint count) {
    AssertLog(sgidx < def()->countSpecs_global());
    pPoolCount[sgidx] += count;
}

////////////////////////////////////////////////////////////////////////////////

uint Raft::getSpecCountByGidx(solver::spec_global_id sgidx) {
    AssertLog(sgidx < def()->countSpecs_global());
    return pPoolCount[sgidx];
}

////////////////////////////////////////////////////////////////////////////////

std::map<steps::triangle_global_id, uint> Raft::getTriSpecCounts(solver::spec_global_id sgidx,
                                                                 const rng::RNGptr rng) {
    std::map<triangle_global_id, uint> tri_counts;
    AssertLog(sgidx < def()->countSpecs_global());

    uint count = pPoolCount[sgidx];
    if (count == 0) {
        return tri_counts;
    }
    uint ntris = pTris_overlap_vec.size();

    // Integer division will give us the floor, which is what we need
    uint count_eq = count / ntris;

    uint count_acc = 0;
    // First assign the equal share
    for (auto const& tri: pTris_overlap_vec) {
        tri_counts[tri] = count_eq;
        count_acc += count_eq;
    }

    // Just in case we hit an equal share by chance
    if (count_acc == count) {
        return tri_counts;
    }

    AssertLog(count_acc < count);

    while (count_acc < count) {
        auto tri = pTris_overlap_vec[rng->get() % ntris];
        tri_counts[tri] += 1;
        count_acc += 1;
    }

    return tri_counts;
}


////////////////////////////////////////////////////////////////////////////////

bool Raft::applyEndoAndDis(double raft_dt) {
    AssertLog(pAppliedEndo == false);

    for (auto& endo: raftendos()) {
        double rate = endo->rate();
        if (rate > 0.0) {
            double endo_dt = patch()->rng()->getExp(rate);
            if (endo_dt < raft_dt) {
                endo->apply();
                // it's possible endo fails if no space for vesicle. If it
                // was succesful it has set this flag
                if (pAppliedEndo) {
                    return true;
                }
            }
        }
    }

    for (auto& dis: raftdiss()) {
        double rate = dis->rate();
        if (rate > 0.0) {
            double dis_dt = patch()->rng()->getExp(rate);
            if (dis_dt < raft_dt) {
                dis->apply();
                return true;
            }
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Raft::setRaftSReacActive(solver::raftsreac_global_id rsridx, bool active) {
    if (active) {
        pRaftSReac_inactive.erase(rsridx);
    } else {
        pRaftSReac_inactive.insert(rsridx);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Raft::getRaftSReacActive(solver::raftsreac_global_id rsridx) const {
    return pRaftSReac_inactive.find(rsridx) == pRaftSReac_inactive.end();
}

////////////////////////////////////////////////////////////////////////////////

void Raft::updImmobility(int mob_upd) {
    if (pImmobility == 0 and mob_upd < 0) {
        std::ostringstream os;
        os << "Negative immobility is not possible for raft. Model error. ";
        ProgErrLog(os.str());
    }

    pImmobility += mob_upd;
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::mpi::tetvesicle
