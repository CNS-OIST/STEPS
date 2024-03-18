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

#include "mpi/tetvesicle/vesreac.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/vessreacdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

static inline double comp_ccst(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

VesReac::VesReac(solver::VesSReacdef* vsrdef, TetRDEF* tet)
    : pVesSReacdef(vsrdef)
    , pTet(tet)
    , pRate_zero(false)
    , pTotal_rate(0.0) {
    AssertLog(pVesSReacdef != nullptr);
    AssertLog(pTet != nullptr);

    pType = KP_VESSREAC;

    // In the situation that the vesicle reaction requires surface, i.e.
    // interacting with species in patches, yet this tet has no surface
    // triangle neighbours then this rate is always 0. Makes sense to flag for
    // this obviously
    double max_distance = def()->max_distance();
    pReqSurface = pVesSReacdef->reqSurface();

    // Kill two birds with one stone here- also set the tri barycentres if we have
    // a max_distance
    bool notri = true;
    for (auto const& tetnexttri: pTet->nexttris()) {
        if (tetnexttri == nullptr) {
            continue;
        }
        if (max_distance > 0.0) {
            pTri_barycentres.push_back(tetnexttri->position());
        }
        notri = false;
    }

    if (pReqSurface && notri) {
        pRate_zero = true;
        // This kproc will never get chosen, no need to set anything else up
        return;
    }
}

////////////////////////////////////////////////////////////////////////////////

VesReac::~VesReac() = default;

////////////////////////////////////////////////////////////////////////////////

void VesReac::checkpoint(std::fstream& cp_file) {
    KProc::checkpoint(cp_file);
    // Note pRate_per_ves and pTotal_rate can simply be set by calling rate()
    // once everything's in place so don't need to be checkpointed
}

////////////////////////////////////////////////////////////////////////////////

void VesReac::restore(std::fstream& cp_file) {
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void VesReac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void VesReac::resetOccupancies() {
    pTet->resetPoolOccupancy();

    for (auto const& tri: pTet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }
        tri->resetPoolOccupancy();
    }
}

////////////////////////////////////////////////////////////////////////////////

void VesReac::setupDeps() {
    if (pRate_zero) {
        return;
    }  // should never get chosen

    KProcPSet updset;

    uint nkprocs = pTet->countKProcs();

    // Search in the tetrahedron.
    for (uint k = 0; k < nkprocs; k++) {
        // Volume species can also affect reacs that belong in this tet
        bool added_k = false;
        for (auto const& s: pVesSReacdef->updColl_O()) {
            if (pTet->KProcDepSpecTet(k, pTet, s)) {
                updset.insert(pTet->getKProc(k));
                added_k = true;
                break;
            }
        }
        if (added_k) {
            continue;
        }
        // Now also vesicle surface species can affect kprocs in tet i.e. vesicle
        // bunding, unbinding
        for (auto const& s: pVesSReacdef->updColl_V()) {
            if (pTet->KProcDepSpecTetVesSurface(k, pTet, s)) {
                updset.insert(pTet->getKProc(k));
                added_k = true;
                break;
            }
        }
        if (added_k) {
            continue;
        }
        // also link species can affect vesicle unbinding
        for (auto const& s: pVesSReacdef->updColl_L()) {
            if (pTet->KProcDepLinkSpecTetVesSurface(k, pTet, s)) {
                updset.insert(pTet->getKProc(k));
                break;
            }
        }
        // Changes on patch triangles can affect other vessreacs
        for (auto const& s: pVesSReacdef->updColl_S()) {
            for (auto const& tri: pTet->nexttris()) {
                if (tri == nullptr) {
                    continue;
                }
                if (pTet->KProcDepSpecTri(k, tri, s)) {
                    updset.insert(pTet->getKProc(k));
                    added_k = true;
                    break;
                }
            }
            if (added_k) {
                break;
            }
        }
    }

    // Volume species changes can also affect surface-type reactions
    // in connected triangles
    for (auto const& tri: pTet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }
        if (pTet->getHost() != tri->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron "
               << pTet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = tri->countKProcs();

        for (uint k = 0; k < nkprocs; k++) {
            bool added_k = false;

            for (auto const& s: pVesSReacdef->updColl_O()) {
                if (tri->KProcDepSpecTet(k, pTet, s)) {
                    updset.insert(tri->getKProc(k));
                    added_k = true;
                    break;
                }
            }

            if (added_k) {
                continue;
            }
            for (auto const& s: pVesSReacdef->updColl_S()) {
                if (tri->KProcDepSpecTri(k, tri, s)) {
                    updset.insert(tri->getKProc(k));
                    break;
                }
            }
        }

        // Tet-owned reactions in other tets that depend on surface specs can also be affected
        TetRDEF* otherTet = (tri->iTet() == pTet) ? tri->oTet() : tri->iTet();
        if (otherTet != nullptr) {
            if (pTet->getHost() != otherTet->getHost()) {
                std::ostringstream os;
                os << "Compartment tetrahedron " << pTet->idx()
                   << " and its patch-separated neighboring tetrahedron " << otherTet->idx()
                   << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            uint nkprocs_ = otherTet->countKProcs();
            for (uint k = 0; k < nkprocs_; k++) {
                for (auto const& s: pVesSReacdef->updColl_S()) {
                    if (otherTet->KProcDepSpecTri(k, tri, s)) {
                        updset.insert(otherTet->getKProc(k));
                        break;
                    }
                }
            }
        }
    }

    updset.insert(this);  // just in case. Shouldn't be necessary now - test

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

double VesReac::rate(TetVesicleRDEF* /*solver*/) {
    if (pRate_zero || inactive()) {
        return 0.0;
    }
    // Kinda tricky but overall rate is summed
    // from 'pairs' where the specific reactants are present in the same
    // tetrahedron thus rate is a total that is summed across all tets that a
    // vesicle currently occupies.

    auto const& vesproxyrefs = pTet->getVesProxyrefs();

    pRate_per_ves.clear();
    pTotal_rate = 0.0;

    // Little time save if total overlap exists.
    // There shouldn't be any species anyway in tet nor any vesicle surface
    if (pTet->vol() == solver::TINY_VOLUME) {
        for (auto const& vesproxy_mapit: vesproxyrefs) {
            auto const& vesproxy_uid = vesproxy_mapit.first;
            AssertLog(pRate_per_ves.find(vesproxy_uid) == pRate_per_ves.end());
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }
        return pTotal_rate;
    }

    double max_distance = def()->max_distance();

    for (auto const& vesproxy_mapit: vesproxyrefs) {
        auto const& vesproxy_uid = vesproxy_mapit.first;
        auto const& vesproxy = vesproxy_mapit.second;
        AssertLog(pRate_per_ves.find(vesproxy_uid) == pRate_per_ves.end());

        // Don't allow if vesicle has been earmarked for exocytosis at next update
        if (vesproxy->exoApplied().valid()) {
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        // TODO // Need to check somewhere these are all the same types of vesicle
        solver::Vesicledef* vdef = vesproxy->def();
        AssertLog(vdef != nullptr);
        solver::vessreac_local_id lidx = vdef->vessreacG2L(pVesSReacdef->gidx());

        // This vessreac is not necessarily defined for the vesicle in question
        if (lidx.unknown()) {
            continue;
        }

        // First check the trigger 'dependency' species. If they are not present on
        // the vesicle surface then rate is 0.

        // Global indices
        const auto& vdep_vec = vdef->vessreac_vdep(lidx);

        bool gotdep = true;
        for (auto s: vdep_vec.range()) {
            uint vdep = vdep_vec[s];
            if (vdep == 0) {
                continue;
            }

            // We got a dependency species. Need to see if it's on the vesicle
            // surface or not.
            uint cnt = vesproxy->getSpecCount_V(s);
            if (cnt >= vdep) {
                continue;
            } else {
                // We are missing a dependency. Rate is 0 for this ves.
                gotdep = false;
                break;
            }
        }

        if (!gotdep) {
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        // If we got here then any vesicle surface species dependencies are present
        // and correct

        // 1: VISIT VESICLE SURFACE SPECIES

        // Need to precompute some triangle information- need to check the minimum
        // distance for vesicle species, only counting those within the distance

        double h_mu = 1.0;

        // Now, the vesicle can have only a subset of species defined for reactions,
        // diffusions but can transport ALL species between membranes. So vesicle
        // holds global indices but they need to be converted to locals or
        // vice-versa

        const auto& lhs_v_vec = vdef->vessreac_lhs_V(lidx);

        for (auto s: lhs_v_vec.range()) {
            uint lhs = lhs_v_vec[s];
            if (lhs == 0) {
                continue;
            }

            solver::spec_global_id spec_gidx = vdef->specL2G(s);

            // get count in the specific tet
            uint cnt = 0;

            if (max_distance > 0.0) {
                const std::vector<std::pair<solver::pointspec_individual_id, math::position_abs>>&
                    spec_positions = vesproxy->getSurfaceSpecPos_absolute(spec_gidx);

                for (auto const& spec_pos: spec_positions) {
                    bool tris_ok = true;
                    for (auto const& tri_pos: pTri_barycentres) {
                        double distance2 = pow(spec_pos.second[0] - tri_pos[0], 2) +
                                           pow(spec_pos.second[1] - tri_pos[1], 2) +
                                           pow(spec_pos.second[2] - tri_pos[2], 2);
                        if (sqrt(distance2) > max_distance) {
                            tris_ok = false;
                        }
                    }
                    if (tris_ok) {
                        cnt += 1;
                    }
                }
            }

            else {
                cnt = vesproxy->getSpecCount_V(spec_gidx);
            }

            if (lhs > cnt) {
                // This is usually a return 0.0 statement, which obviously doesn't make
                // sense now since we're in a loop over vesicles
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
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        // Visit any link species

        const auto& lhs_l_vec = vdef->vessreac_lhs_L(lidx);

        // Now, the vesicle can have only a subset of species defined for reactions,
        // diffusions but can transport ALL species between membranes. So vesicle
        // holds global indices but they need to be converted to locals or
        // vice-versa

        for (auto l: lhs_l_vec.range()) {
            uint lhs = lhs_l_vec[l];
            if (lhs == 0) {
                continue;
            }

            // get count in the specific tet
            uint cnt = vesproxy->getLinkSpecCount_V(l);

            if (lhs > cnt) {
                // This is usually a return 0.0 statement, which obviously doesn't make
                // sense now since we're in a loop over vesicles
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
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        // 2: OUTER TET VOLUME SPECIES

        // Find species in 'outer' tet.

        // Global indices
        const auto& lhs_o_vec = vdef->vessreac_lhs_O(lidx);

        // Local indices
        const auto& cnt_o_vec = pTet->pools();

        for (auto s: lhs_o_vec.range()) {
            uint lhs = lhs_o_vec[s];
            if (lhs == 0) {
                continue;
            }
            // Need to convert species global index to local for the tets (same for
            // all)
            solver::spec_local_id spec_olidx = pTet->compdef()->specG2L(s);

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
        if (h_mu == 0.0) {
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        // 3: INNER  VOLUME SPECIES
        // We no longer allow on lhs

        // 4: ALL SPECIES FROM CONNECTED TRIANGLES TO TET

        // Problem here is that we have potentially 4 different triangles in which
        // we can have any species we need. Any triangle can add to the h_mu for
        // this tet, but if there are none then h_mu becomes zero. Correct thing to
        // do is first check if there are any lhs for surface,
        // - if so then loop through tris and sum the count, use that count to set
        // h_mu,
        // - repeat for other species.
        // May not be optimal, but the only way I can currently think to do it.

        // All species, global indices
        const auto& lhs_s_vec = vdef->vessreac_lhs_S(lidx);

        for (auto sg: lhs_s_vec.range()) {
            uint lhs = lhs_s_vec[sg];
            if (lhs == 0) {
                continue;
            }
            // Let's compute the sum of h_mus in each tri
            double h_mu_mult = 0;

            for (auto const& tetnexttri: pTet->nexttris()) {
                if (tetnexttri == nullptr) {
                    continue;
                }

                // We need local patch indices
                solver::spec_local_id spec_slidx = (tetnexttri)->patchdef()->specG2L(sg);
                auto cnt = tetnexttri->pools()[spec_slidx];
                double tri_h_mu = 1.0;

                if (lhs > cnt) {
                    continue;
                }
                switch (lhs) {
                case 4: {
                    tri_h_mu *= static_cast<double>(cnt - 3);
                }
                case 3: {
                    tri_h_mu *= static_cast<double>(cnt - 2);
                }
                case 2: {
                    tri_h_mu *= static_cast<double>(cnt - 1);
                }
                case 1: {
                    tri_h_mu *= static_cast<double>(cnt);
                    break;
                }
                default: {
                    AssertLog(0);
                    return 0.0;
                }
                }
                h_mu_mult += tri_h_mu;
            }

            h_mu *= h_mu_mult;
        }
        if (h_mu == 0.0) {
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        double ccst = comp_ccst(kcst(), pTet->vol(), pVesSReacdef->order());

        // Multiply with scaled reaction constant.
        double rate = h_mu * ccst;
        pRate_per_ves[vesproxy_uid] = rate;
        pTotal_rate += rate;
    }

    return pTotal_rate;
}

////////////////////////////////////////////////////////////////////////////////

void VesReac::apply(const rng::RNGptr& rng,
                    double /*dt*/,
                    double /*simtime*/,
                    double period,
                    TetVesicleRDEF* /*solver*/) {
    AssertLog(!pRate_zero);

    double selector = rng->getUnfII() * pTotal_rate;

    // Select by rate
    solver::vesicle_individual_id vesproxy_uid;
    double accum = 0.0;
    for (auto const& vesr: pRate_per_ves) {
        vesproxy_uid = vesr.first;
        accum += vesr.second;
        if (selector < accum) {
            break;
        }
    }

    AssertLog(vesproxy_uid.valid());

    VesProxy* vesproxy = pTet->getVesProxyref(vesproxy_uid);

    solver::Vesicledef* vdef = vesproxy->def();
    solver::vessreac_local_id lidx = vdef->vessreacG2L(pVesSReacdef->gidx());

    // If we have defined a rate for this vesicle then the vessreac better belong
    // to that vesicle
    AssertLog(lidx.valid());

    // Update vesicle counts. Do this first because it can cause the apply to fail
    std::map<solver::spec_global_id, int> spec_upd;
    const auto& upd_v_vec = vdef->vessreac_upd_V(lidx);

    for (auto v: upd_v_vec.range()) {
        int upd = upd_v_vec[v];
        if (upd == 0) {
            continue;
        }
        // In terms of reactions, the vesicle deals with local indices, but
        // can transport other species. Need to convert to GLOBAL species indices
        solver::spec_global_id spec_gidx = vdef->specL2G(v);

        spec_upd[spec_gidx] = upd;
    }

    bool added = vesproxy->addSurfSpecs(spec_upd);

    if (!added) {
        CLOG(WARNING, "general_log") << "Failed to apply vesicle surface reaction. Couldn't find "
                                        "position for vesicle surface spec\n";
        return;
    }


    // Update link species

    // Here, if we have an update, we should have a -1 update and a +1 update.
    // This gives us the indices to chance in the linkspecies.
    const auto& upd_l_vec = vdef->vessreac_upd_L(lidx);

    solver::linkspec_global_id posupd_specgidx;
    solver::linkspec_global_id negupd_specgidx;

    for (auto l: upd_l_vec.range()) {
        int upd = upd_l_vec[l];
        if (upd == 0) {
            continue;
        }
        // l is global index
        if (upd == 1) {
            posupd_specgidx = l;
        } else if (upd == -1) {
            negupd_specgidx = l;
        } else {
            std::ostringstream os;
            os << "Unexpected update of link species of >1 or <-1. ";
            ProgErrLog(os.str());
        }
    }

    if ((posupd_specgidx.valid() && negupd_specgidx.unknown()) ||
        (posupd_specgidx.unknown() && negupd_specgidx.valid())) {
        std::ostringstream os;
        os << "Unexpected update of link species. Link species cannot be created "
              "or destroyed: only change identity. ";
        ProgErrLog(os.str());
    }

    if (posupd_specgidx.valid()) {
        vesproxy->changeLinkSpecGidx(negupd_specgidx, posupd_specgidx);
    }

    // Update outer tet pools.
    const auto& upd_o_vec = vdef->vessreac_upd_O(lidx);
    const auto& cnt_o_vec = pTet->pools();
    for (auto s: upd_o_vec.range()) {
        // Best to do this check first because spec may be undefined in tet if this is 0
        int upd = upd_o_vec[s];
        if (upd == 0) {
            continue;
        }
        // s is global index. Need local index for tet
        solver::spec_local_id spec_tetlidx = pTet->compdef()->specG2L(s);

        AssertLog(spec_tetlidx.valid());

        if (pTet->clamped(spec_tetlidx) == true) {
            continue;
        }

        int nc = static_cast<int>(cnt_o_vec[spec_tetlidx]) + upd;
        AssertLog(nc >= 0);
        pTet->setCount(spec_tetlidx, static_cast<uint>(nc), period);
    }

    // Update inner vescile pools.
    const auto& upd_i_vec = vdef->vessreac_upd_I(lidx);
    for (auto s: upd_i_vec.range()) {
        uint upd = upd_i_vec[s];
        if (upd == 0) {
            continue;
        }
        vesproxy->incSpecCount_I(s, upd);
    }

    // Update patch pools.
    // If this tet's rate is non-zero then its triangle neighbours should be
    // populated correctly. Have to select one from potentially 4 for the update.

    const auto& upd_s_vec = vdef->vessreac_upd_S(lidx);
    // Need to go spec by spec, selecting all tris that can be updated, then
    // randomly choosing one
    for (auto s: upd_s_vec.range()) {
        int upd = upd_s_vec[s];
        if (upd == 0) {
            continue;
        }
        // Need some of this stuff outside the loop

        std::map<uint, uint> update_tris;

        for (uint i = 0; i <= 3; ++i) {
            TriRDEF* tri = pTet->nextTri(i);
            if (tri == nullptr) {
                continue;
            }
            // S is global index. Need local index for tri
            solver::spec_local_id spec_trilidx = tri->patchdef()->specG2L(s);

            // If we got here upd is non-zero so spec should be defined
            AssertLog(spec_trilidx.valid());

            uint cnt_s = tri->pools()[spec_trilidx];

            // If we're adding species they can be added to any triangle
            if (upd > 0) {
                // Even if clamped it's important to add this for selection, even though
                // it's unchanging
                if (tri->clamped(spec_trilidx)) {
                    update_tris[i] = cnt_s;
                } else {
                    update_tris[i] = cnt_s + upd;
                }
            } else {
                // Removing (a) molecule(s)

                // Only add this tri if there are enough molecules
                if (cnt_s >= static_cast<uint>(-upd)) {
                    // Even if clamped it's important to add this for selection, even
                    // though it's unchanging
                    if (tri->clamped(spec_trilidx)) {
                        update_tris[i] = cnt_s;
                    } else {
                        update_tris[i] = cnt_s + upd;
                    }
                }
            }
        }

        uint tupdate_size = update_tris.size();

        // Should have gotten at least one triangle if rate is non-zero, or only +ve
        // update
        AssertLog(tupdate_size > 0);

        // Now select from vector and update
        uint rIndex = rng->get() % tupdate_size;

        uint idx = 0;
        // Just for testing/debugging
        bool set = false;
        for (auto const& upd_it: update_tris) {
            if (idx == rIndex) {
                TriRDEF* tri = pTet->nextTri(upd_it.first);
                // Have to get this again here because tris may be in different patches
                solver::spec_local_id spec_trilidx = tri->patchdef()->specG2L(s);
                tri->setCount(spec_trilidx, upd_it.second, period);
                set = true;

                break;
            }

            idx += 1;
        }
        // Just for testing/debugging
        AssertLog(set);
    }

    // Finally set immobilization or not
    vesproxy->updImmobility(pVesSReacdef->immobility());

    rExtent++;

    // Extent held also at the def level since these objects come in and out of
    // existence
    pVesSReacdef->incExtent();
}

}  // namespace steps::mpi::tetvesicle
