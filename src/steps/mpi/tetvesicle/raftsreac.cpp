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

#include "mpi/tetvesicle/raftsreac.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/raftsreacdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

static inline double comp_ccst(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_area(double kcst, double area, uint order) {
    double ascale = area * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    double ccst = kcst * pow(ascale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac::RaftSReac(solver::RaftSReacdef* rsrdef, TriRDEF* tri)
    : pRaftSReacdef(rsrdef)
    , pTri(tri)
    , pTotal_rate(0.0) {
    AssertLog(pRaftSReacdef != nullptr);
    AssertLog(pRaftSReacdef != nullptr);

    pType = KP_RAFTSREAC;
}

////////////////////////////////////////////////////////////////////////////////

RaftSReac::~RaftSReac() = default;

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::checkpoint(std::fstream& cp_file) {
    KProc::checkpoint(cp_file);
    // Note pRate_per_raft and pTotal_rate can simply be set by calling rate()
    // once everything's in place so don't need to be checkpointed
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::restore(std::fstream& cp_file) {
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

// This function now calculates the deps on the fly depending which tet was
// applied to ALTERNATIVELY could pre-calculate deps for every tetrahedron and
// store in a map?
void RaftSReac::setupDeps() {
    // The 'update columns' contain anything changed by the reaction,
    // whether created or destroyed.
    // Anything that changes on the raft surface can affect other raft reactions
    // BUT all raft stuff gets a kind of global update after every SSA step at the
    // moment Anything that changes in the volume can affect any kproc in that
    // volume or surface-reaction types in a connected triangle

    KProcPSet updset;

    TetRDEF* itet = pTri->iTet();
    TetRDEF* otet = pTri->oTet();

    uint nkprocs = pTri->countKProcs();
    for (uint sk = 0; sk < nkprocs; sk++) {
        for (auto const& spec: pRaftSReacdef->updColl_S()) {
            if (pTri->KProcDepSpecTri(sk, pTri, spec)) {
                updset.insert(pTri->getKProc(sk));
                break;  // exit the species loop because we can only add this kp once
            }
        }
    }

    // Do this in a separate loop. Alternative to set a flag if kp has already
    // been added for pRaftSReacdef->updColl_S to avoid going through
    // pRaftSReacdef->updColl_Rs too when unneccesary
    for (uint sk = 0; sk < nkprocs; sk++) {
        for (auto const& spec: pRaftSReacdef->updColl_Rs()) {
            if (pTri->KProcDepSpecTriRaftSurface(sk, pTri, spec)) {
                updset.insert(pTri->getKProc(sk));
                break;  // exit the species loop because we can only add this kp once
            }
        }
    }

    if (itet != nullptr) {
        if (pTri->getHost() != itet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
            os << itet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = itet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            bool added_k = false;
            for (auto const& spec: pRaftSReacdef->updColl_I()) {
                if (itet->KProcDepSpecTet(k, itet, spec)) {
                    updset.insert(itet->getKProc(k));
                    added_k = true;
                    break;  // exit the species loop because we can only add this kp once
                }
            }
            if (added_k) {
                continue;
            }
            // Changes in surface species affect reactions like vessreacs
            for (auto const& spec: pRaftSReacdef->updColl_S()) {
                if (itet->KProcDepSpecTri(k, pTri, spec)) {
                    updset.insert(itet->getKProc(k));
                    break;  // exit the species loop because we can only add this kp once
                }
            }
        }

        for (auto const& tri_i: itet->nexttris()) {
            if (tri_i == nullptr) {
                continue;
            }

            if (tri_i->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << tri_i->idx() << " and its compartment tetrahedron ";
                os << itet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = tri_i->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                for (auto const& spec: pRaftSReacdef->updColl_I()) {
                    if (tri_i->KProcDepSpecTet(sk, itet, spec)) {
                        updset.insert(tri_i->getKProc(sk));
                        break;  // exit the species loop because we can only add this kp once
                    }
                }
            }
        }
    }

    if (otet != nullptr) {
        if (pTri->getHost() != otet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
            os << otet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = otet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            bool added_k = false;
            for (auto const& spec: pRaftSReacdef->updColl_O()) {
                if (otet->KProcDepSpecTet(k, otet, spec)) {
                    updset.insert(otet->getKProc(k));
                    added_k = true;
                    break;  // exit the species loop because we can only add this kp once
                }
            }
            if (added_k) {
                continue;
            }
            // Changes in surface species affect reactions like vessreacs
            for (auto const& spec: pRaftSReacdef->updColl_S()) {
                if (otet->KProcDepSpecTri(k, pTri, spec)) {
                    updset.insert(otet->getKProc(k));
                    break;  // exit the species loop because we can only add this kp once
                }
            }
        }

        for (auto const& tri_o: otet->nexttris()) {
            if (tri_o == nullptr) {
                continue;
            }
            if (tri_o->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << tri_o->idx() << " and its compartment tetrahedron ";
                os << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = tri_o->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                for (auto const& spec: pRaftSReacdef->updColl_O()) {
                    if (tri_o->KProcDepSpecTet(sk, otet, spec)) {
                        updset.insert(tri_o->getKProc(sk));
                        break;  // exit the species loop because we can only add this kp once
                    }
                }
            }
        }
    }

    updset.insert(this);  // just in case. Shouldn't be necessary now - test

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

// TODO, do I Need to add a new function here like depSpecRaft which returns
// true of this raft surface reaction contains spec gidx on lhs in raft surface

////////////////////////////////////////////////////////////////////////////////

double RaftSReac::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // overall rate is summed
    // from 'pairs' where the specific reactants are present in the same
    // tetrahedron thus rate is a total that is summed across all tets that a raft
    // currently occupies.

    // First get the current list of rafts overlapping the parent tri
    auto raftproxyrefs = pTri->getRaftProxyrefs();

    pRate_per_raft.clear();
    pTotal_rate = 0.0;

    for (auto const& raftproxyref: raftproxyrefs) {
        auto const& raftproxy_uid = raftproxyref.first;
        auto const& raftproxy = raftproxyref.second;
        solver::Raftdef* rdef = raftproxy->def();
        solver::raftsreac_local_id lidx = rdef->raftsreacG2L(pRaftSReacdef->gidx());

        // This raftsreac is not necessarily defined for the raft in question
        if (lidx.unknown()) {
            continue;
        }

        if (!raftproxy->getRaftSReacActive(pRaftSReacdef->gidx())) {
            continue;
        }

        // Check dependency species first

        // Global indices
        const auto& rsdep_vec = rdef->raftsreac_rsdep(lidx);

        bool gotdep = true;
        for (auto s: rsdep_vec.range()) {
            uint rsdep = rsdep_vec[s];
            if (rsdep == 0) {
                continue;
            }
            // We got a dependency species. Need to check raft surface
            uint cnt = raftproxy->getSpecCountByGidx(s);
            if (cnt >= rsdep) {
                continue;
            } else {
                // We are missing a dependency. Rate is 0 for this raft.
                gotdep = false;
                break;
            }
        }

        if (!gotdep) {
            pRate_per_raft[raftproxy_uid] = 0.0;
            continue;
        }

        const auto& anti_rsdep_vec = rdef->raftsreac_anti_rsdep(lidx);

        bool gotantidep = false;

        for (auto s: anti_rsdep_vec.range()) {
            uint anti_rsdep = anti_rsdep_vec[s];
            if (anti_rsdep == 0) {
                continue;
            }
            // We got an anti-dependency species. Need to check raft surface
            // Note: no check here s is actually global index
            uint cnt = raftproxy->getSpecCountByGidx(s);
            if (cnt < anti_rsdep) {
                continue;
            } else {
                // We have hit an anti-dependency. Rate is 0 for this raft.
                gotantidep = true;
                break;
            }
        }

        if (gotantidep) {
            pRate_per_raft[raftproxy_uid] = 0.0;
            continue;
        }

        // 1: VISIT RAFT SURFACE SPECIES

        double h_mu = 1.0;

        // Now, the raft can have only a subset of species defined for reactions,
        // diffusions but can transport ALL species between membranes. So raft holds
        // global indices but they need to be converted to locals or vice-versa

        // We have to do a boolean test for lhs because rdef->countSpecs contains
        // both lhs and rhs species, and we only want to know lhs for scaling (see
        // note below this loop)
        const auto& lhs_rs_vec = rdef->raftsreac_lhs_Rs(lidx);

        for (auto s: lhs_rs_vec.range()) {
            uint lhs = lhs_rs_vec[s];
            if (lhs == 0) {
                continue;
            }

            solver::spec_global_id spec_gidx = rdef->specL2G(s);

            // get count in the specific tet
            uint cnt = raftproxy->pools_global()[spec_gidx.get()];

            if (lhs > cnt) {
                h_mu = 0.0;
                break;
            }
            switch (lhs) {
            case 4: {
                h_mu *= static_cast<double>(cnt - 3);
            }
            case 3: {
                h_mu *= static_cast<double>(cnt - 2);
            }
            case 2: {
                h_mu *= static_cast<double>(cnt - 1);
            }
            case 1: {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default: {
                AssertLog(0);
                return 0.0;
            }
            }
        }

        if (h_mu == 0.0) {
            pRate_per_raft[raftproxy_uid] = 0.0;
            continue;
        }

        const auto& lhs_s_vec = rdef->raftsreac_lhs_S(lidx);
        const auto& cnt_s_vec = pTri->pools();

        for (auto s: lhs_s_vec.range()) {
            uint lhs = lhs_s_vec[s];
            if (lhs == 0) {
                continue;
            }

            // The following relies on s being global index- how to check this is the case??
            // Need to convert species global index to local for the tris
            solver::spec_local_id spec_slidx = pTri->patchdef()->specG2L(s);

            uint cnt = cnt_s_vec[spec_slidx];

            if (lhs > cnt) {
                h_mu = 0.0;
                break;
            }
            switch (lhs) {
            case 4: {
                h_mu *= static_cast<double>(cnt - 3);
            }
            case 3: {
                h_mu *= static_cast<double>(cnt - 2);
            }
            case 2: {
                h_mu *= static_cast<double>(cnt - 1);
            }
            case 1: {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default: {
                AssertLog(0);
                return 0.0;
            }
            }
        }

        if (h_mu == 0.0) {
            pRate_per_raft[raftproxy_uid] = 0.0;
            continue;
        }

        // We only need to visit tets if it's not surf-surf
        if (pRaftSReacdef->surf_surf() == false) {
            if (def()->outside()) {
                // 2: OUTER TET VOLUME SPECIES

                // Find species in 'outer' tet.

                // Global indices
                const auto& lhs_o_vec = rdef->raftsreac_lhs_O(lidx);

                // Local indices
                const auto& cnt_o_vec = pTri->oTet()->pools();

                for (auto s: lhs_o_vec.range()) {
                    uint lhs = lhs_o_vec[s];
                    if (lhs == 0) {
                        continue;
                    }
                    // Need to convert species global index to local for the tets (same
                    // for all). How to check this is indeed global index?
                    solver::spec_local_id spec_olidx = pTri->oTet()->compdef()->specG2L(s);

                    uint cnt = cnt_o_vec[spec_olidx];

                    if (lhs > cnt) {
                        h_mu = 0.0;
                        break;
                    }
                    switch (lhs) {
                    case 4: {
                        h_mu *= static_cast<double>(cnt - 3);
                    }
                    case 3: {
                        h_mu *= static_cast<double>(cnt - 2);
                    }
                    case 2: {
                        h_mu *= static_cast<double>(cnt - 1);
                    }
                    case 1: {
                        h_mu *= static_cast<double>(cnt);
                        break;
                    }
                    default: {
                        AssertLog(0);
                        return 0.0;
                    }
                    }
                }
            }

            if (def()->inside()) {
                // 3: INNER  VOLUME SPECIES

                // Global indices
                const auto& lhs_i_vec = rdef->raftsreac_lhs_I(lidx);

                // Local indices
                const auto& cnt_i_vec = pTri->iTet()->pools();

                for (auto s: lhs_i_vec.range()) {
                    uint lhs = lhs_i_vec[s];
                    if (lhs == 0) {
                        continue;
                    }
                    // Need to convert species global index to local for the tets (same
                    // for all). How to check s is global??
                    solver::spec_local_id spec_ilidx = pTri->iTet()->compdef()->specG2L(s);

                    uint cnt = cnt_i_vec[spec_ilidx];

                    if (lhs > cnt) {
                        h_mu = 0.0;
                        break;
                    }

                    switch (lhs) {
                    case 4: {
                        h_mu *= static_cast<double>(cnt - 3);
                    }
                    case 3: {
                        h_mu *= static_cast<double>(cnt - 2);
                    }
                    case 2: {
                        h_mu *= static_cast<double>(cnt - 1);
                    }
                    case 1: {
                        h_mu *= static_cast<double>(cnt);
                        break;
                    }
                    default: {
                        AssertLog(0);
                        return 0.0;
                    }
                    }
                }
            }
        }

        if (h_mu == 0.0) {
            pRate_per_raft[raftproxy_uid] = 0.0;
            continue;
        }

        double ccst = 0.0;

        if (pRaftSReacdef->surf_surf() == false) {
            double vol;
            if (pRaftSReacdef->inside() == true) {
                AssertLog(pTri->iTet() != nullptr);
                vol = pTri->iTet()->vol();
            } else {
                AssertLog(pTri->oTet() != nullptr);
                vol = pTri->oTet()->vol();
            }

            ccst = comp_ccst(kcst(), vol, pRaftSReacdef->order());
        } else {
            double area;
            area = pTri->area();
            ccst = comp_ccst_area(kcst(), area, pRaftSReacdef->order());
        }

        // Multiply with scaled reaction constant.
        double rate = h_mu * ccst;
        pRate_per_raft[raftproxy_uid] = rate;
        pTotal_rate += rate;
    }

    return pTotal_rate;
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::apply(const rng::RNGptr& rng,
                      double /*dt*/,
                      double /*simtime*/,
                      double period,
                      TetVesicleRDEF* /*solver*/) {
    // Pick a tet for the application based on pools, make the changes, then
    // calculate deps based on that tet (using function getDeps) and return that
    // update vector

    double selector = rng->getUnfII() * pTotal_rate;

    // Select by rate
    solver::raft_individual_id raftproxy_uid;
    double accum = 0.0;
    for (auto const& rft: pRate_per_raft) {
        raftproxy_uid = rft.first;
        accum += rft.second;
        if (selector < accum) {
            break;
        }
    }

    AssertLog(raftproxy_uid.valid());

    RaftProxy* raftproxy = pTri->getRaftProxyref(raftproxy_uid);

    solver::Raftdef* rdef = raftproxy->def();
    solver::raftsreac_local_id lidx = rdef->raftsreacG2L(pRaftSReacdef->gidx());

    // Update raft counts.
    std::map<uint, int> spec_upd;

    const auto& upd_rs_vec = rdef->raftsreac_upd_Rs(lidx);
    const auto& cnt_rs_vec = raftproxy->pools_global();

    for (auto s: upd_rs_vec.range()) {
        int upd = upd_rs_vec[s];
        if (upd == 0) {
            continue;
        }
        solver::spec_global_id spec_gidx = rdef->specL2G(s);
        int nc = static_cast<int>(cnt_rs_vec[spec_gidx.get()]) + upd;
        AssertLog(nc >= 0);
        raftproxy->setSpecCountByLidx(s, static_cast<uint>(nc));
    }

    // Update tri species
    // This is global indices again
    const auto& upd_s_vec = rdef->raftsreac_upd_S(lidx);
    const auto& cnt_s_vec = pTri->pools();
    for (auto s: upd_s_vec.range()) {
        // upd vec comes from raft and so is global indices
        int upd = upd_s_vec[s];
        if (upd == 0) {
            continue;
        }
        // S is global index. Need local index for tri. Should exist if upd is non-zero
        solver::spec_local_id spec_trilidx = pTri->patchdef()->specG2L(s);
        AssertLog(spec_trilidx.valid());

        if (pTri->clamped(spec_trilidx) == true) {
            continue;
        }

        int nc = static_cast<int>(cnt_s_vec[spec_trilidx]) + upd;
        AssertLog(nc >= 0);
        pTri->setCount(spec_trilidx, static_cast<uint>(nc), period);
    }

    // Update outer tet pools.
    TetRDEF* otet = pTri->oTet();
    if (otet != nullptr) {
        const auto& upd_o_vec = rdef->raftsreac_upd_O(lidx);
        const auto& cnt_o_vec = otet->pools();
        for (auto s: upd_o_vec.range()) {
            // Best to do this check first because spec may be undefined in tet if this is 0
            int upd = upd_o_vec[s];
            if (upd == 0) {
                continue;
            }
            // S is global index. Need local index for tet
            solver::spec_local_id spec_tetlidx = otet->compdef()->specG2L(s);

            AssertLog(spec_tetlidx.valid());

            if (otet->clamped(spec_tetlidx) == true) {
                continue;
            }

            int nc = static_cast<int>(cnt_o_vec[spec_tetlidx]) + upd;
            AssertLog(nc >= 0);
            otet->setCount(spec_tetlidx, static_cast<uint>(nc), period);
        }
    }

    // Update inner tet pools.
    TetRDEF* itet = pTri->iTet();
    if (itet != nullptr) {
        const auto& upd_i_vec = rdef->raftsreac_upd_I(lidx);
        const auto& cnt_i_vec = itet->pools();
        for (auto s: upd_i_vec.range()) {
            int upd = upd_i_vec[s];
            if (upd == 0) {
                continue;
            }

            // S is global index. Need local index for tet
            solver::spec_local_id spec_tetlidx = itet->compdef()->specG2L(s);

            AssertLog(spec_tetlidx.valid());

            if (itet->clamped(spec_tetlidx) == true) {
                continue;
            }

            int nc = static_cast<int>(cnt_i_vec[spec_tetlidx]) + upd;
            AssertLog(nc >= 0);
            itet->setCount(spec_tetlidx, static_cast<uint>(nc), period);
        }
    }

    // Finally set immobilization or not
    raftproxy->updImmobility(pRaftSReacdef->immobility());

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

void RaftSReac::resetOccupancies() {
    pTri->resetPoolOccupancy();

    // Update inner tet pools.
    TetRDEF* itet = pTri->iTet();
    if (itet != nullptr) {
        itet->resetPoolOccupancy();
    }

    TetRDEF* otet = pTri->oTet();
    if (otet != nullptr) {
        otet->resetPoolOccupancy();
    }
}

}  // namespace steps::mpi::tetvesicle
