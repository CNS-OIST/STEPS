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
#include "sreac.hpp"
#include "math/constants.hpp"
#include "tetexact.hpp"
#include "tri.hpp"
#include "wmvol.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

#include "util/checkpointing.hpp"

namespace steps::tetexact {


////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_vol(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // I.H 5/1/2011 Removed this strange special behaviour for zero-order
    // if (o1 < 0) o1 = 0;

    // BUGFIX, IH 16/8/2010. Incorrect for zero-order reactions:
    // rate was the same for all tets, not a fraction of the total volume
    // return kcst * pow(vscale, static_cast<double>(-o1));

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    // Special case for zero-order reactions right now: - removed, not necessary
    // now units are correct
    // if (order == 0) ccst *= (area/patcharea);

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

SReac::SReac(solver::SReacdef* srdef, Tri* tri)
    : pSReacdef(srdef)
    , pTri(tri)
    , pCcst(0.0)
    , pKcst(0.0) {
    AssertLog(pSReacdef != nullptr);
    AssertLog(pTri != nullptr);

    solver::sreac_local_id lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
    double kcst = pTri->patchdef()->kcst(lsridx);
    pKcst = kcst;

    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside()) {
            AssertLog(pTri->iTet() != nullptr);
            vol = pTri->iTet()->vol();
        } else {
            AssertLog(pTri->oTet() != nullptr);
            vol = pTri->oTet()->vol();
        }

        pCcst = comp_ccst_vol(kcst, vol, pSReacdef->order());
    } else {
        double area;
        area = pTri->area();
        pCcst = comp_ccst_area(kcst, area, pSReacdef->order());
    }

    AssertLog(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

SReac::~SReac() = default;

////////////////////////////////////////////////////////////////////////////////

void SReac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    _resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::_resetCcst() {
    solver::sreac_local_id lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
    double kcst = pTri->patchdef()->kcst(lsridx);
    pKcst = kcst;

    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside()) {
            AssertLog(pTri->iTet() != nullptr);
            vol = pTri->iTet()->vol();
        } else {
            AssertLog(pTri->oTet() != nullptr);
            vol = pTri->oTet()->vol();
        }

        pCcst = comp_ccst_vol(kcst, vol, pSReacdef->order());
    } else {
        double area;
        area = pTri->area();
        pCcst = comp_ccst_area(kcst, area, pSReacdef->order());
    }

    AssertLog(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;

    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside()) {
            AssertLog(pTri->iTet() != nullptr);
            vol = pTri->iTet()->vol();
        } else {
            AssertLog(pTri->oTet() != nullptr);
            vol = pTri->oTet()->vol();
        }

        pCcst = comp_ccst_vol(k, vol, pSReacdef->order());
    } else {
        double area;
        area = pTri->area();
        pCcst = comp_ccst_area(k, area, pSReacdef->order());
    }

    AssertLog(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setupDeps() {
    // For all non-zero entries gidx in SReacDef's UPD_S:
    //   Perform depSpecTri(gidx,tri()) for:
    //     All kproc's of tri()
    //     NOTE: we currently don't test kproc's of inner and external tet,
    //     because that's not strictly necessary.
    //
    // If inner tetrahedron exists:
    //   For all non-zero entries gidx in SReacDef's UPD_I:
    //     Perform depSpecTet(gidx,itet) for:
    //       All kproc's of itet
    //       All kproc's of triangles next to itet
    //
    // If outer tetrahedron exists:
    //   Similar to inner tet.
    //
    // All dependencies are first collected into a std::set, to sort them
    // and to eliminate duplicates. At the end of the routine, they are
    // copied into the vector that will be returned during execution.

    WmVol* itet = pTri->iTet();
    WmVol* otet = pTri->oTet();

    KProcPSet updset;
    for (auto const& k: pTri->kprocs()) {
        for (auto const& spec: pSReacdef->updColl_S()) {
            if (k->depSpecTri(spec, pTri)) {
                updset.insert(k);
            }
        }
    }

    if (itet != nullptr) {
        for (auto const& k: itet->kprocs()) {
            for (auto const& spec: pSReacdef->updColl_I()) {
                if (k->depSpecTet(spec, itet)) {
                    updset.insert(k);
                }
            }
        }


        for (auto const& tri: itet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }

            for (auto const& k: tri->kprocs()) {
                for (auto const& spec: pSReacdef->updColl_I()) {
                    if (k->depSpecTet(spec, itet)) {
                        updset.insert(k);
                    }
                }
            }
        }
    }

    if (otet != nullptr) {
        for (auto const& k: otet->kprocs()) {
            for (auto const& spec: pSReacdef->updColl_O()) {
                if (k->depSpecTet(spec, otet)) {
                    updset.insert(k);
                }
            }
        }

        for (auto const& tri: otet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }

            for (auto const& k: tri->kprocs()) {
                for (auto const& spec: pSReacdef->updColl_O()) {
                    if (k->depSpecTet(spec, otet)) {
                        updset.insert(k);
                    }
                }
            }
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecTet(solver::spec_global_id gidx, WmVol* tet) {
    // We need to check whether the tet is inside or outside.
    //   -> If inside: check dependency using SReacDef's I_DEP
    //   -> If outside: check dependency using SReacDef's O_DEP
    //   -> If neither, return.
    //
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp

    if (tet == pTri->iTet()) {
        return (pSReacdef->dep_I(gidx) != solver::DEP_NONE);
    } else if (tet == pTri->oTet()) {
        return (pSReacdef->dep_O(gidx) != solver::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecTri(solver::spec_global_id gidx, Tri* triangle) {
    if (triangle != pTri) {
        return false;
    }
    return (pSReacdef->dep_S(gidx) != solver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

double SReac::rate(Tetexact* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // First we compute the combinatorial part.
    //   1/ for the surface part of the stoichiometry
    //   2/ for the inner or outer volume part of the stoichiometry,
    //      depending on whether the sreac is inner() or outer()
    // Then we multiply with mesoscopic constant.

    solver::Patchdef* pdef = pTri->patchdef();
    solver::sreac_local_id lidx = pdef->sreacG2L(pSReacdef->gidx());

    double h_mu = 1.0;

    const auto& lhs_s_vec = pdef->sreac_lhs_S(lidx);
    const auto& cnt_s_vec = pTri->pools();
    for (auto s: lhs_s_vec.range()) {
        uint lhs = lhs_s_vec[s];
        if (lhs == 0) {
            continue;
        }
        uint cnt = cnt_s_vec[s];
        if (lhs > cnt) {
            return 0.0;
        }
        switch (lhs) {
        case 4: {
            h_mu *= static_cast<double>(cnt - 3);
        }
            STEPS_FALLTHROUGH;
        case 3: {
            h_mu *= static_cast<double>(cnt - 2);
        }
            STEPS_FALLTHROUGH;
        case 2: {
            h_mu *= static_cast<double>(cnt - 1);
        }
            STEPS_FALLTHROUGH;
        case 1: {
            h_mu *= static_cast<double>(cnt);
            break;
        }
        default: {
            AssertLog(0);
        }
        }
    }

    if (pSReacdef->inside()) {
        const auto& lhs_i_vec = pdef->sreac_lhs_I(lidx);
        auto const& cnt_i_vec = pTri->iTet()->pools();
        for (auto s: lhs_i_vec.range()) {
            uint lhs = lhs_i_vec[s];
            if (lhs == 0) {
                continue;
            }
            uint cnt = cnt_i_vec[s];
            if (lhs > cnt) {
                return 0.0;
            }
            switch (lhs) {
            case 4: {
                h_mu *= static_cast<double>(cnt - 3);
            }
                STEPS_FALLTHROUGH;
            case 3: {
                h_mu *= static_cast<double>(cnt - 2);
            }
                STEPS_FALLTHROUGH;
            case 2: {
                h_mu *= static_cast<double>(cnt - 1);
            }
                STEPS_FALLTHROUGH;
            case 1: {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default: {
                AssertLog(false);
            }
            }
        }
    } else if (pSReacdef->outside()) {
        const auto& lhs_o_vec = pdef->sreac_lhs_O(lidx);
        auto const& cnt_o_vec = pTri->oTet()->pools();
        for (auto s: lhs_o_vec.range()) {
            uint lhs = lhs_o_vec[s];
            if (lhs == 0) {
                continue;
            }
            uint cnt = cnt_o_vec[s];
            if (lhs > cnt) {
                return 0.0;
            }
            switch (lhs) {
            case 4: {
                h_mu *= static_cast<double>(cnt - 3);
            }
                STEPS_FALLTHROUGH;
            case 3: {
                h_mu *= static_cast<double>(cnt - 2);
            }
                STEPS_FALLTHROUGH;
            case 2: {
                h_mu *= static_cast<double>(cnt - 1);
            }
                STEPS_FALLTHROUGH;
            case 1: {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default: {
                AssertLog(0);
            }
            }
        }
    }

    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& SReac::apply(const rng::RNGptr& /*rng*/, double dt, double simtime) {
    solver::Patchdef* pdef = pTri->patchdef();
    solver::sreac_local_id lidx = pdef->sreacG2L(pSReacdef->gidx());

    // Update triangle pools.
    const auto& upd_s_vec = pdef->sreac_upd_S(lidx);
    const auto& cnt_s_vec = pTri->pools();

    // First tell the triangles of any channel states relating to ohmic currents
    // that have been changed. If a change has occurred then the triangle will
    // calculate the current that passed through the
    // channel since the last change.
    // NOTE: This must clearly come BEFORE the actual update happens
    for (auto oc: solver::ohmiccurr_local_id::range(pdef->countOhmicCurrs())) {
        // Patchdef returns local index
        solver::spec_local_id cs_lidx = pdef->ohmiccurr_chanstate(oc);
        if (pTri->clamped(cs_lidx)) {
            continue;
        }
        int upd = upd_s_vec[cs_lidx];
        if (upd == 0) {
            continue;
        }
        // Now here a channel state related to an ohmic current has changed
        // it's number: tell the triangle
        pTri->setOCchange(oc, cs_lidx, dt, simtime);
    }

    for (auto s: upd_s_vec.range()) {
        if (pTri->clamped(s)) {
            continue;
        }
        int upd = upd_s_vec[s];
        if (upd == 0) {
            continue;
        }
        int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        AssertLog(nc >= 0);
        pTri->setCount(s, static_cast<uint>(nc));
    }

    // Update inner tet pools.
    WmVol* itet = pTri->iTet();
    if (itet != nullptr) {
        const auto& upd_i_vec = pdef->sreac_upd_I(lidx);
        auto const& cnt_i_vec = itet->pools();
        for (auto s: upd_i_vec.range()) {
            if (itet->clamped(s)) {
                continue;
            }
            int upd = upd_i_vec[s];
            if (upd == 0) {
                continue;
            }
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            AssertLog(nc >= 0);
            itet->setCount(s, static_cast<uint>(nc));
        }
    }

    // Update outer tet pools.
    WmVol* otet = pTri->oTet();
    if (otet != nullptr) {
        const auto& upd_o_vec = pdef->sreac_upd_O(lidx);
        auto const& cnt_o_vec = otet->pools();
        for (auto s: upd_o_vec.range()) {
            if (otet->clamped(s)) {
                continue;
            }
            int upd = upd_o_vec[s];
            if (upd == 0) {
                continue;
            }
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            AssertLog(nc >= 0);
            otet->setCount(s, static_cast<uint>(nc));
        }
    }

    rExtent++;

    return pUpdVec;
}

}  // namespace steps::tetexact
