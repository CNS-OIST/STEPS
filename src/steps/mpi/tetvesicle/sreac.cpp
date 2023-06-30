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

#include "mpi/tetvesicle/sreac.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

// logging
#include "easylogging++.h"

namespace steps::mpi::tetvesicle {

static inline double comp_ccst_vol(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    auto o1 = static_cast<int>(order) - 1;

    // I.H 5/1/2011 Removed this strange special behaviour for zero-order
    // if (o1 < 0) o1 = 0;

    // BUGFIX, IH 16/8/2010. Incorrect for zero-order reactions:
    // rate was the same for all tets, not a fraction of the total volume
    // return kcst * pow(vscale, static_cast<double>(-o1));

    auto ccst = kcst * pow(vscale, static_cast<double>(-o1));

    // Special case for zero-order reactions right now: - removed, not necessary
    // now units are correct
    // if (order == 0) ccst *= (area/patcharea);

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_area(double kcst, double area, uint order) {
    double ascale = area * math::AVOGADRO;
    auto o1 = static_cast<int>(order) - 1;

    double ccst = kcst * pow(ascale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

SReac::SReac(solver::SReacdef* srdef, TriRDEF* tri)
    : pSReacdef(srdef)
    , pTri(tri)
    , pCcst(0.0)
    , pKcst(0.0) {
    AssertLog(pSReacdef != nullptr);
    AssertLog(pTri != nullptr);

    pType = KP_SREAC;

    solver::sreac_local_id lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
    double kcst = pTri->patchdef()->kcst(lsridx);
    pKcst = kcst;

    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside() == true) {
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
    resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::resetCcst() {
    solver::sreac_local_id lsridx = pTri->patchdef()->sreacG2L(pSReacdef->gidx());
    double kcst = pTri->patchdef()->kcst(lsridx);
    pKcst = kcst;

    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside() == true) {
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
        if (pSReacdef->inside() == true) {
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

    AssertLog(pTri->getInHost());

    TetRDEF* itet = pTri->iTet();
    TetRDEF* otet = pTri->oTet();

    KProcPSet updset;

    uint nkprocs = pTri->countKProcs();
    // check if sk KProc in pTri depends on spec in pTri
    for (uint sk = 0; sk < nkprocs; sk++) {
        for (auto const& spec: pSReacdef->updColl_S()) {
            if (pTri->KProcDepSpecTri(sk, pTri, spec)) {
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

        // Checking Kprocs in inner tet that might depend on spec change in that tet
        nkprocs = itet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            bool added_k = false;
            for (auto const& spec: pSReacdef->updColl_I()) {
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
            for (auto const& spec: pSReacdef->updColl_S()) {
                if (itet->KProcDepSpecTri(k, pTri, spec)) {
                    updset.insert(itet->getKProc(k));
                    break;  // exit the species loop because we can only add this kp once
                }
            }
        }

        // Checking Kprocs in triangles attached to inner tet that might
        // depend on spec change in that tet (such as sreacs)
        for (auto const& tri: itet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }
            if (tri->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron ";
                os << itet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = tri->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                for (auto const& spec: pSReacdef->updColl_I()) {
                    if (tri->KProcDepSpecTet(sk, itet, spec)) {
                        updset.insert(tri->getKProc(sk));
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

        // Checking Kprocs in outer tet that might depend on spec change in that tet
        nkprocs = otet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            bool added_k = false;
            for (auto const& spec: pSReacdef->updColl_O()) {
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
            for (auto const& spec: pSReacdef->updColl_S()) {
                if (otet->KProcDepSpecTri(k, pTri, spec)) {
                    updset.insert(otet->getKProc(k));
                    break;  // exit the species loop because we can only add this kp once
                }
            }
        }

        // Checking Kprocs in triangles attached to outer tet that might
        // depend on spec change in that tet (such as sreacs)
        for (auto const& tri: otet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }
            if (tri->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron ";
                os << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = tri->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                for (auto const& spec: pSReacdef->updColl_O()) {
                    if (tri->KProcDepSpecTet(sk, otet, spec)) {
                        updset.insert(tri->getKProc(sk));
                        break;  // exit the species loop because we can only add this kp once
                    }
                }
            }
        }
    }

    localUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

double SReac::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // Have to do this business because volumes (and areas) might
    // have changed by vesicle occupancy
    if (pSReacdef->surf_surf() == false) {
        double vol;
        if (pSReacdef->inside() == true) {
            AssertLog(pTri->iTet() != nullptr);
            vol = pTri->iTet()->vol();
        } else {
            AssertLog(pTri->oTet() != nullptr);
            vol = pTri->oTet()->vol();
        }

        pCcst = comp_ccst_vol(pKcst, vol, pSReacdef->order());
    } else {
        double area;
        area = pTri->area();
        pCcst = comp_ccst_area(pKcst, area, pSReacdef->order());
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
            AssertLog(false);
        }
        }
    }

    if (pSReacdef->inside()) {
        const auto& lhs_i_vec = pdef->sreac_lhs_I(lidx);
        const auto& cnt_i_vec = pTri->iTet()->pools();
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
                AssertLog(false);
            }
            }
        }
    } else if (pSReacdef->outside()) {
        const auto& lhs_o_vec = pdef->sreac_lhs_O(lidx);
        const auto& cnt_o_vec = pTri->oTet()->pools();
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

    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::apply(const rng::RNGptr& /*rng*/,
                  double dt,
                  double simtime,
                  double period,
                  TetVesicleRDEF* /*solver*/) {
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
    uint nocs = pdef->countOhmicCurrs();
    for (auto oc: solver::ohmiccurr_local_id::range(nocs)) {
        // Patchdef returns local index
        solver::spec_local_id cs_lidx = pdef->ohmiccurr_chanstate(oc);
        if (pTri->clamped(cs_lidx) == true) {
            continue;
        }
        int upd = upd_s_vec[cs_lidx];
        if (upd == 0) {
            continue;
        }  // Now here a channel state related to an ohmic current has changed
        // it's number: tell the triangle
        pTri->setOCchange(oc, cs_lidx, dt, simtime);
    }

    for (auto s: upd_s_vec.range()) {
        if (pTri->clamped(s) == true) {
            continue;
        }
        int upd = upd_s_vec[s];
        if (upd == 0) {
            continue;
        }
        int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        AssertLog(nc >= 0);
        pTri->setCount(s, static_cast<uint>(nc), period);
    }

    // Update inner tet pools.
    TetRDEF* itet = pTri->iTet();
    if (itet != nullptr) {
        const auto& upd_i_vec = pdef->sreac_upd_I(lidx);
        const auto& cnt_i_vec = itet->pools();
        for (auto s: upd_i_vec.range()) {
            if (itet->clamped(s) == true) {
                continue;
            }
            int upd = upd_i_vec[s];
            if (upd == 0) {
                continue;
            }
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            AssertLog(nc >= 0);
            itet->setCount(s, static_cast<uint>(nc), period);
        }
    }

    // Update outer tet pools.
    TetRDEF* otet = pTri->oTet();
    if (otet != nullptr) {
        const auto& upd_o_vec = pdef->sreac_upd_O(lidx);
        const auto& cnt_o_vec = otet->pools();
        for (auto s: upd_o_vec.range()) {
            if (otet->clamped(s) == true) {
                continue;
            }
            int upd = upd_o_vec[s];
            if (upd == 0) {
                continue;
            }
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            AssertLog(nc >= 0);
            otet->setCount(s, static_cast<uint>(nc), period);
        }
    }

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& SReac::getRemoteUpdVec(int /*direction*/) const {
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& SReac::getLocalUpdVec(int /*direction*/) const {
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void SReac::resetOccupancies() {
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
