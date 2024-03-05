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
#include "comp.hpp"
#include "math/constants.hpp"
#include "wmrssa.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::wmrssa {


////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order) {
    double vscale = 1.0e3 * vol * math::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(solver::Reacdef* rdef, Comp* comp)
    : pReacdef(rdef)
    , pComp(comp)
    , pCcst(0.0) {
    assert(pReacdef != nullptr);
    assert(pComp != nullptr);
    solver::reac_local_id lridx = pComp->def()->reacG2L(pReacdef->gidx());
    double kcst = pComp->def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
    assert(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac() = default;

////////////////////////////////////////////////////////////////////////////////

void Reac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    util::checkpoint(cp_file, pPropensityLB);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    util::restore(cp_file, pPropensityLB);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::active() const {
    solver::reac_local_id lridx = pComp->def()->reacG2L(defr()->gidx());
    return pComp->def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::reset() {
    resetExtent();
    solver::reac_local_id lridx = pComp->def()->reacG2L(defr()->gidx());
    pComp->def()->setActive(lridx, true);
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

void Reac::resetCcst() {
    solver::reac_local_id lridx = pComp->def()->reacG2L(pReacdef->gidx());
    double kcst = pComp->def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp->def()->vol(), pReacdef->order());
    assert(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setupDeps() {
    SchedIDXSet updset;

    // Search in local compartment.
    for (auto const& k: pComp->kprocs()) {
        for (auto const& s: defr()->updColl()) {
            if (k->depSpecComp(s, pComp)) {
                updset.insert(k->schedIDX());
            }
        }
    }

    // Search in neighbouring patches.
    for (auto const& p: pComp->ipatches()) {
        for (auto const& k: p->kprocs()) {
            for (auto const& s: defr()->updColl()) {
                if (k->depSpecComp(s, pComp)) {
                    updset.insert(k->schedIDX());
                }
            }
        }
    }
    for (auto const& p: pComp->opatches()) {
        for (auto const& k: p->kprocs()) {
            for (auto const& s: defr()->updColl()) {
                if (k->depSpecComp(s, pComp)) {
                    updset.insert(k->schedIDX());
                }
            }
        }
    }

    schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecComp(solver::spec_global_id gidx, Comp* comp) {
    if (pComp != comp) {
        return false;
    }
    return defr()->dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecPatch(solver::spec_global_id /*gidx*/, Patch* /*patch*/) {
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double Reac::rate(wmrssa::PropensityRSSA prssa) {
    if (inactive()) {
        return 0.0;
    }

    if (prssa == wmrssa::BOUNDS) {
        pPropensityLB = rate(wmrssa::LOWERBOUND);
    }

    // Prefetch some variables.
    solver::Compdef* cdef = pComp->def();
    const auto& lhs_vec = cdef->reac_lhs(cdef->reacG2L(defr()->gidx()));
    const auto& cnt_vec = pComp->pools(prssa);

    // Compute combinatorial part.
    double h_mu = 1.0;
    for (auto pool: lhs_vec.range()) {
        uint lhs = lhs_vec[pool];
        if (lhs == 0) {
            continue;
        }
        auto cnt = static_cast<uint>(cnt_vec[pool]);
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
            AssertLog(false);
        }
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& Reac::apply() {
    SchedIDXSet updset;
    // bool returnUpdVec = false;
    solver::Compdef* cdef = pComp->def();
    const auto& local = cdef->pools();
    solver::reac_local_id l_ridx = cdef->reacG2L(defr()->gidx());
    const auto& upd_vec = cdef->reac_upd(l_ridx);
    for (auto i: upd_vec.range()) {
        if (cdef->clamped(i)) {
            continue;
        }
        int j = upd_vec[i];
        if (j == 0) {
            continue;
        }
        int nc = static_cast<int>(local[i]) + j;
        cdef->setCount(i, static_cast<double>(nc));
        if (pComp->isOutOfBound(i, nc)) {
            for (const auto& dependentReac: pComp->getSpecUpdKProcs(i)) {
                updset.insert(dependentReac->schedIDX());
            }
            // returnUpdVec = true;
        }
    }
    rExtent++;
    schedIDXSet_To_Vec(updset, pUpdVec);
    return pUpdVec;  // returnUpdVec ? pUpdVec : emptyVec;
}

}  // namespace steps::wmrssa
