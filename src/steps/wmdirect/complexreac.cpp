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
#include "wmdirect/complexreac.hpp"
#include "math/constants.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"
#include "util/error.hpp"
#include "wmdirect/comp.hpp"
#include "wmdirect/kproc.hpp"
#include "wmdirect/wmdirect.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::wmdirect {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order) {
    const double vscale = 1.0e3 * vol * math::AVOGADRO;
    const int o1 = static_cast<int>(order) - 1;

    // IMPORTANT: Now treating zero-order reaction units correctly, i.e. as
    // M/s not /s
    // if (o1 < 0) o1 = 0;

    return kcst * std::pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

ComplexReac::ComplexReac(solver::ComplexReacdef& rdef, Comp& comp)
    : pComplexReacdef(rdef)
    , pComp(comp) {
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

ComplexReac::~ComplexReac() = default;

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexReac::active() const {
    solver::complexreac_local_id lridx = pComp.def()->complexreacG2L(defcr().gidx());
    return pComp.def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::reset() {
    resetExtent();
    solver::complexreac_local_id lridx = pComp.def()->complexreacG2L(defcr().gidx());
    pComp.def()->setComplexReacActive(lridx, true);
    resetCcst();
    for (auto& v: candidates) {
        v.second.reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::resetCcst() {
    solver::complexreac_local_id lridx = pComp.def()->complexreacG2L(pComplexReacdef.gidx());
    const double kcst = pComp.def()->kcst(lridx);
    pCcst = comp_ccst(kcst, pComp.def()->vol(), pComplexReacdef.order());
    AssertLog(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexReac::setupDeps() {
    SchedIDXSet updset;

    // Search in local compartment.
    for (auto const& k: pComp.kprocs()) {
        for (auto const& s: defcr().updColl()) {
            if (k->depSpecComp(s, &pComp)) {
                updset.insert(k->schedIDX());
            }
        }
        for (const auto& v: defcr().complexUPDMAP()) {
            for (auto sus: v.second) {
                if (k->depComplexComp(v.first, sus, &pComp)) {
                    updset.insert(k->schedIDX());
                }
            }
        }
    }

    // Search in neighbouring patches.
    for (auto const& p: pComp.ipatches()) {
        for (auto const& k: p->kprocs()) {
            for (auto const& s: defcr().updColl()) {
                if (k->depSpecComp(s, &pComp)) {
                    updset.insert(k->schedIDX());
                }
            }
            for (const auto& v: defcr().complexUPDMAP()) {
                for (auto sus: v.second) {
                    if (k->depComplexComp(v.first, sus, &pComp)) {
                        updset.insert(k->schedIDX());
                    }
                }
            }
        }
    }
    for (auto const& p: pComp.opatches()) {
        for (auto const& k: p->kprocs()) {
            for (auto const& s: defcr().updColl()) {
                if (k->depSpecComp(s, &pComp)) {
                    updset.insert(k->schedIDX());
                }
            }
            for (const auto& v: defcr().complexUPDMAP()) {
                for (auto sus: v.second) {
                    if (k->depComplexComp(v.first, sus, &pComp)) {
                        updset.insert(k->schedIDX());
                    }
                }
            }
        }
    }

    schedIDXSet_To_Vec(updset, pUpdVec);

    // Setup complex events
    for (auto ue: defcr().updEvents()) {
        const solver::complex_global_id cmplxId(ue->complexIdx());
        auto it = candidates.find(cmplxId);
        if (it == candidates.end()) {
            it = candidates.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
        }
        it->second.addEvent(ue, *pComp.def());
    }
    for (auto de: defcr().delEvents()) {
        const solver::complex_global_id cmplxId(de->complexIdx());
        auto it = candidates.find(cmplxId);
        if (it == candidates.end()) {
            it = candidates.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
        }
        it->second.addEvent(de, *pComp.def());
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexReac::depSpecComp(solver::spec_global_id gidx, Comp* comp) {
    return &pComp == comp and defcr().dep(gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexReac::depComplexComp(solver::complex_global_id gidx,
                                 solver::complex_substate_id sus,
                                 Comp* comp) {
    return &pComp == comp and defcr().complexdep(gidx, sus);
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexReac::depSpecPatch(solver::spec_global_id /*gidx*/, Patch* /*patch*/) {
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double ComplexReac::rate() const {
    if (inactive()) {
        return 0.0;
    }

    solver::Compdef* cdef = pComp.def();

    // Get the rates for the complex reaction part
    double cmult = 1.0;
    for (auto& v: candidates) {
        cmult *= v.second.rateMult(cdef->complexStates(v.first));
    }

    if (cmult == 0.0) {
        return 0.0;
    }

    // Prefetch some variables.
    const auto& lhs_vec = cdef->complexreac_lhs(cdef->complexreacG2L(defcr().gidx()));
    const auto& cnt_vec = cdef->pools();

    // Compute species combinatorial part.
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

    // Multiply with scaled reaction constant.
    return h_mu * cmult * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& ComplexReac::apply() {
    // Species part
    solver::Compdef* cdef = pComp.def();
    const auto& local = cdef->pools();
    const auto& upd_vec = cdef->complexreac_upd(cdef->complexreacG2L(defcr().gidx()));

    for (auto i: upd_vec.range()) {
        if (cdef->clamped(i)) {
            continue;
        }
        const int j = upd_vec[i];
        if (j == 0) {
            continue;
        }
        const int nc = static_cast<int>(local[i]) + j;
        cdef->setCount(i, static_cast<double>(nc));
    }

    // Complex part
    auto& rng = defcr().statedef().rng();
    // Updates
    for (auto& v: candidates) {
        for (auto& event: v.second.selectEvents(rng)) {
            if (event.first->type() == solver::UPDEvent) {
                const auto ev = std::dynamic_pointer_cast<const solver::ComplexUpdateEventdef>(
                    event.first);
                const auto& state = cdef->complexStates(v.first).at(event.second);
                cdef->updateComplex(v.first, event.second, ev->getUpdate(state, rng));
            } else {
                // Deletions
                cdef->removeComplex(v.first, event.second);
            }
        }
    }
    // Creations
    for (auto& ce: defcr().creEvents()) {
        cdef->addComplex(ce->complexIdx(), pComp.solver()->getNextComplexInd(), ce->init());
    }

    rExtent++;
    return pUpdVec;
}

}  // namespace steps::wmdirect
