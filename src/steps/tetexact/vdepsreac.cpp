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
#include "vdepsreac.hpp"
#include "math/constants.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "tri.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetexact {


////////////////////////////////////////////////////////////////////////////////

VDepSReac::VDepSReac(solver::VDepSReacdef* vdsrdef, Tri* tri)
    : pVDepSReacdef(vdsrdef)
    , pTri(tri)
    , pScaleFactor(0.0) {
    AssertLog(pVDepSReacdef != nullptr);
    AssertLog(pTri != nullptr);

    if (pVDepSReacdef->surf_surf() == false) {
        double vol;
        if (pVDepSReacdef->inside()) {
            AssertLog(pTri->iTet() != nullptr);
            vol = pTri->iTet()->vol();
        } else {
            AssertLog(pTri->oTet() != nullptr);
            vol = pTri->oTet()->vol();
        }
        double vscale = 1.0e3 * vol * math::AVOGADRO;
        uint order = pVDepSReacdef->order();
        double o1 = static_cast<double>(order) - 1.0;

        pScaleFactor = pow(vscale, -o1);
    } else {
        double ascale = pTri->area() * math::AVOGADRO;
        uint order = pVDepSReacdef->order();
        double o1 = static_cast<double>(order) - 1.0;

        pScaleFactor = pow(ascale, -o1);
    }

    AssertLog(pScaleFactor > 0.0);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac::~VDepSReac() = default;

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::checkpoint(std::fstream& cp_file) {
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::restore(std::fstream& cp_file) {
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void VDepSReac::setupDeps() {
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
        for (auto const& spec: pVDepSReacdef->updcoll_S()) {
            if (k->depSpecTri(spec, pTri)) {
                updset.insert(k);
            }
        }
    }

    if (itet != nullptr) {
        for (auto const& k: itet->kprocs()) {
            for (auto const& spec: pVDepSReacdef->updcoll_I()) {
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
                for (auto const& spec: pVDepSReacdef->updcoll_I()) {
                    if (k->depSpecTet(spec, itet)) {
                        updset.insert(k);
                    }
                }
            }
        }
    }

    if (otet != nullptr) {
        for (auto const& k: otet->kprocs()) {
            for (auto const& spec: pVDepSReacdef->updcoll_O()) {
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
                for (auto const& spec: pVDepSReacdef->updcoll_O()) {
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

bool VDepSReac::depSpecTet(solver::spec_global_id gidx, WmVol* tet) {
    // We need to check whether the tet is inside or outside.
    //   -> If inside: check dependency using VDepSReacDef's I_DEP
    //   -> If outside: check dependency using VDepSReacDef's O_DEP
    //   -> If neither, return.
    //
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp

    if (tet == pTri->iTet()) {
        return pVDepSReacdef->dep_I(gidx) != solver::DEP_NONE;
    } else if (tet == pTri->oTet()) {
        return pVDepSReacdef->dep_O(gidx) != solver::DEP_NONE;
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool VDepSReac::depSpecTri(solver::spec_global_id gidx, Tri* triangle) {
    if (triangle != pTri) {
        return false;
    }
    return pVDepSReacdef->dep_S(gidx) != solver::DEP_NONE;
}

////////////////////////////////////////////////////////////////////////////////

double VDepSReac::rate(Tetexact* solver) {
    if (inactive()) {
        return 0.0;
    }

    // First we compute the combinatorial part.
    //   1/ for the surface part of the stoichiometry
    //   2/ for the inner or outer volume part of the stoichiometry,
    //      depending on whether the sreac is inner() or outer()
    // Then we multiply with mesoscopic constant.

    solver::Patchdef* pdef = pTri->patchdef();
    solver::vdepsreac_local_id lidx = pdef->vdepsreacG2L(pVDepSReacdef->gidx());

    double h_mu = 1.0;

    const auto& lhs_s_vec = pdef->vdepsreac_lhs_S(lidx);
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
            AssertLog(false);
        }
        }
    }

    if (pVDepSReacdef->inside()) {
        const auto& lhs_i_vec = pdef->vdepsreac_lhs_I(lidx);
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
    } else if (pVDepSReacdef->outside()) {
        const auto& lhs_o_vec = pdef->vdepsreac_lhs_O(lidx);
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
                AssertLog(false);
            }
            }
        }
    }

    double v = solver->getTriV(pTri->idx());
    double k = pVDepSReacdef->getVDepK(v);

    return h_mu * k * pScaleFactor;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& VDepSReac::apply(const rng::RNGptr& /*rng*/, double dt, double simtime) {
    // NOTE: simtime is BEFORE the update has taken place

    solver::Patchdef* pdef = pTri->patchdef();
    solver::vdepsreac_local_id lidx = pdef->vdepsreacG2L(pVDepSReacdef->gidx());

    const auto& upd_s_vec = pdef->vdepsreac_upd_S(lidx);
    const auto& cnt_s_vec = pTri->pools();

    // First tell the triangles of any channel states relating to ohmic currents
    // that have been changed. If a change has occured then the triangle will
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

    // Update triangle pools.
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
        const auto& upd_i_vec = pdef->vdepsreac_upd_I(lidx);
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
        const auto& upd_o_vec = pdef->vdepsreac_upd_O(lidx);
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
