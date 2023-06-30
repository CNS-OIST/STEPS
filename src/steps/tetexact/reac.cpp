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
#include "reac.hpp"

#include "math/constants.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "wmvol.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetexact {


////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order, double /*compvol*/) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    double ccst = kcst * pow(vscale, static_cast<double>(-o1));

    return ccst;
}

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(solver::Reacdef* rdef, WmVol* tet)
    : pReacdef(rdef)
    , pTet(tet)
    , pCcst(0.0)
    , pKcst(0.0) {
    AssertLog(pReacdef != nullptr);
    AssertLog(pTet != nullptr);

    solver::reac_local_id lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void Reac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pKcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pKcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    _resetCcst();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::_resetCcst() {
    solver::reac_local_id lridx = pTet->compdef()->reacG2L(pReacdef->gidx());
    double kcst = pTet->compdef()->kcst(lridx);
    // Also reset kcst
    pKcst = kcst;
    pCcst = comp_ccst(kcst, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setKcst(double k) {
    AssertLog(k >= 0.0);
    pKcst = k;
    pCcst = comp_ccst(k, pTet->vol(), pReacdef->order(), pTet->compdef()->vol());
    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setupDeps() {
    KProcPSet updset;

    // Search in local tetrahedron.
    for (auto const& k: pTet->kprocs()) {
        for (auto const& s: pReacdef->UPD_Coll()) {
            if (k->depSpecTet(s, pTet)) {
                // updset.insert((*k)->getSSARef());
                updset.insert(k);
            }
        }
    }

    for (auto const& tri: pTet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }

        for (auto const& k: tri->kprocs()) {
            for (auto const& s: pReacdef->UPD_Coll()) {
                if (k->depSpecTet(s, pTet) == true) {
                    // updset.insert((*k)->getSSARef());
                    updset.insert(k);
                }
            }
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());
    // pUpdObjVec.assign(updset_obj.begin(), updset_obj.end());
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecTet(solver::spec_global_id gidx, WmVol* tet) {
    if (pTet != tet) {
        return false;
    }
    return pReacdef->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecTri(solver::spec_global_id /*gidx*/, Tri* /*tri*/) {
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double Reac::rate(Tetexact* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // Prefetch some variables.
    solver::Compdef* cdef = pTet->compdef();
    const auto& lhs_vec = cdef->reac_lhs(cdef->reacG2L(pReacdef->gidx()));
    auto const& cnt_vec = pTet->pools();

    // Compute combinatorial part.
    double h_mu = 1.0;
    for (auto pool: lhs_vec.range()) {
        uint lhs = lhs_vec[pool];
        if (lhs == 0) {
            continue;
        }
        uint cnt = cnt_vec[pool];
        if (lhs > cnt) {
            h_mu = 0.0;
            break;
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
            return 0.0;
        }
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& Reac::apply(const rng::RNGptr& /*rng*/,
                                       double /*dt*/,
                                       double /*simtime*/) {
    auto const& local = pTet->pools();
    solver::Compdef* cdef = pTet->compdef();
    solver::reac_local_id l_ridx = cdef->reacG2L(pReacdef->gidx());
    const auto& upd_vec = cdef->reac_upd(l_ridx);
    for (auto i: upd_vec.range()) {
        if (pTet->clamped(i)) {
            continue;
        }
        int j = upd_vec[i];
        if (j == 0) {
            continue;
        }
        int nc = static_cast<int>(local[i]) + j;
        pTet->setCount(i, static_cast<uint>(nc));
    }
    rExtent++;
    return pUpdVec;
}

}  // namespace steps::tetexact
