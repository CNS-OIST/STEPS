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

#include "wmdirect/complexsreac.hpp"

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

static inline double comp_ccst_vol(double kcst, double vol, uint order) {
    const double vscale = 1.0e3 * vol * math::AVOGADRO;
    const int o1 = static_cast<int>(order) - 1;
    return kcst * std::pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst_area(double kcst, double area, uint order) {
    const double ascale = area * math::AVOGADRO;
    const int o1 = static_cast<int>(order) - 1;
    return kcst * std::pow(ascale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

ComplexSReac::ComplexSReac(solver::ComplexSReacdef& rdef, Patch& patch)
    : pComplexSReacdef(rdef)
    , pPatch(patch) {
    resetCcst();
}

////////////////////////////////////////////////////////////////////////////////

ComplexSReac::~ComplexSReac() = default;

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pCcst);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::restore(std::fstream& cp_file) {
    util::restore(cp_file, pCcst);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReac::active() const {
    solver::complexsreac_local_id lridx = pPatch.def()->complexsreacG2L(defcsr().gidx());
    return pPatch.def()->active(lridx);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::reset() {
    resetExtent();
    solver::complexsreac_local_id lridx = pPatch.def()->complexsreacG2L(defcsr().gidx());
    pPatch.def()->setComplexSReacActive(lridx, true);
    resetCcst();
    for (auto& cand: candidates) {
        for (auto& v: cand.second) {
            v.second.reset();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::resetCcst() {
    solver::complexsreac_local_id lridx = pPatch.def()->complexsreacG2L(pComplexSReacdef.gidx());
    const double kcst = pPatch.def()->kcst(lridx);

    if (not defcsr().surf_surf()) {
        double vol;
        if (defcsr().inside()) {
            AssertLog(pPatch.iComp() != nullptr);
            vol = pPatch.iComp()->def()->vol();
        } else {
            AssertLog(pPatch.oComp() != nullptr);
            vol = pPatch.oComp()->def()->vol();
        }
        pCcst = comp_ccst_vol(kcst, vol, defcsr().order());
    } else {
        pCcst = comp_ccst_area(kcst, pPatch.def()->area(), defcsr().order());
    }

    AssertLog(pCcst >= 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void ComplexSReac::setupDeps() {
    Comp* icomp = pPatch.iComp();
    Comp* ocomp = pPatch.oComp();
    SchedIDXSet updset;

    // Search in local patch.
    for (auto const& k: pPatch.kprocs()) {
        for (auto const& s: defcsr().updColl_S()) {
            if (k->depSpecPatch(s, &pPatch)) {
                updset.insert(k->schedIDX());
            }
        }
        for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_SURF)) {
            for (auto sus: v.second) {
                if (k->depComplexPatch(v.first, sus, &pPatch)) {
                    updset.insert(k->schedIDX());
                }
            }
        }
    }

    if (icomp != nullptr) {
        for (auto const& k: icomp->kprocs()) {
            for (auto const& spec: defcsr().updColl_I()) {
                if (k->depSpecComp(spec, icomp)) {
                    updset.insert(k->schedIDX());
                }
            }
            for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_IN)) {
                for (auto sus: v.second) {
                    if (k->depComplexComp(v.first, sus, icomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
            }
        }
        for (auto const& ip: icomp->ipatches()) {
            for (auto const& k: ip->kprocs()) {
                for (auto const& spec: defcsr().updColl_I()) {
                    if (k->depSpecComp(spec, icomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
                for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_IN)) {
                    for (auto sus: v.second) {
                        if (k->depComplexComp(v.first, sus, icomp)) {
                            updset.insert(k->schedIDX());
                        }
                    }
                }
            }
        }
        for (auto const& op: icomp->opatches()) {
            for (auto const& k: op->kprocs()) {
                for (auto const& spec: defcsr().updColl_I()) {
                    if (k->depSpecComp(spec, icomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
                for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_IN)) {
                    for (auto sus: v.second) {
                        if (k->depComplexComp(v.first, sus, icomp)) {
                            updset.insert(k->schedIDX());
                        }
                    }
                }
            }
        }
    }

    if (ocomp != nullptr) {
        for (auto const& k: ocomp->kprocs()) {
            for (auto const& spec: defcsr().updColl_O()) {
                if (k->depSpecComp(spec, ocomp)) {
                    updset.insert(k->schedIDX());
                }
            }
            for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_OUT)) {
                for (auto sus: v.second) {
                    if (k->depComplexComp(v.first, sus, ocomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
            }
        }
        for (auto const& ip: ocomp->ipatches()) {
            for (auto const& k: ip->kprocs()) {
                for (auto const& spec: defcsr().updColl_O()) {
                    if (k->depSpecComp(spec, ocomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
                for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_OUT)) {
                    for (auto sus: v.second) {
                        if (k->depComplexComp(v.first, sus, ocomp)) {
                            updset.insert(k->schedIDX());
                        }
                    }
                }
            }
        }
        for (auto const& op: ocomp->opatches()) {
            for (auto const& k: op->kprocs()) {
                for (auto const& spec: defcsr().updColl_O()) {
                    if (k->depSpecComp(spec, ocomp)) {
                        updset.insert(k->schedIDX());
                    }
                }
                for (const auto& v: defcsr().complexUPDMAP(model::ComplexLocation::PATCH_OUT)) {
                    for (auto sus: v.second) {
                        if (k->depComplexComp(v.first, sus, ocomp)) {
                            updset.insert(k->schedIDX());
                        }
                    }
                }
            }
        }
    }

    schedIDXSet_To_Vec(updset, pUpdVec);

    // Setup complex events
    {
        const model::ComplexLocation loc = model::ComplexLocation::PATCH_SURF;
        auto& cands = candidates[loc];
        for (auto& ue: defcsr().updEvents(loc)) {
            solver::complex_global_id cmplxId(ue->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(ue, *pPatch.def());
        }
        for (auto& de: defcsr().delEvents(loc)) {
            solver::complex_global_id cmplxId(de->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(de, *pPatch.def());
        }
    }
    if (icomp != nullptr) {
        const model::ComplexLocation loc = model::ComplexLocation::PATCH_IN;
        auto& cands = candidates[loc];
        for (auto& ue: defcsr().updEvents(loc)) {
            solver::complex_global_id cmplxId(ue->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(ue, *icomp->def());
        }
        for (auto& de: defcsr().delEvents(loc)) {
            solver::complex_global_id cmplxId(de->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(de, *icomp->def());
        }
    }
    if (ocomp != nullptr) {
        const model::ComplexLocation loc = model::ComplexLocation::PATCH_OUT;
        auto& cands = candidates[loc];
        for (auto& ue: defcsr().updEvents(loc)) {
            solver::complex_global_id cmplxId(ue->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(ue, *ocomp->def());
        }
        for (auto& de: defcsr().delEvents(loc)) {
            solver::complex_global_id cmplxId(de->complexIdx());
            auto it = cands.find(cmplxId);
            if (it == cands.end()) {
                it = cands.emplace(cmplxId, ComplexLHSCandidates(cmplxId)).first;
            }
            it->second.addEvent(de, *ocomp->def());
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReac::depSpecComp(solver::spec_global_id gidx, Comp* comp) {
    if (comp == pPatch.iComp()) {
        return (defcsr().dep_I(gidx) != solver::DEP_NONE);
    } else if (comp == pPatch.oComp()) {
        return (defcsr().dep_O(gidx) != solver::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReac::depComplexComp(solver::complex_global_id gidx,
                                  solver::complex_substate_id sus,
                                  Comp* comp) {
    if (comp == pPatch.iComp()) {
        return defcsr().complexdep(model::ComplexLocation::PATCH_IN, gidx, sus);
    } else if (comp == pPatch.oComp()) {
        return defcsr().complexdep(model::ComplexLocation::PATCH_OUT, gidx, sus);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReac::depSpecPatch(solver::spec_global_id gidx, Patch* patch) {
    return patch == &pPatch and defcsr().dep_S(gidx) != solver::DEP_NONE;
}

////////////////////////////////////////////////////////////////////////////////

bool ComplexSReac::depComplexPatch(solver::complex_global_id gidx,
                                   solver::complex_substate_id sus,
                                   Patch* patch) {
    return patch == &pPatch and defcsr().complexdep(model::ComplexLocation::PATCH_SURF, gidx, sus);
}

////////////////////////////////////////////////////////////////////////////////

double ComplexSReac::rate() const {
    if (inactive()) {
        return 0.0;
    }

    const solver::Patchdef* pdef = pPatch.def();
    solver::complexsreac_local_id lidx = pdef->complexsreacG2L(defcsr().gidx());

    double h_mu = 1.0;
    double cmult = 1.0;

    // Surface part
    const auto& lhs_s_vec = pdef->complexsreac_lhs_S(lidx);
    const auto& cnt_s_vec = pdef->pools();
    for (auto s: lhs_s_vec.range()) {
        const uint lhs = lhs_s_vec[s];
        if (lhs == 0) {
            continue;
        }
        const auto cnt = static_cast<uint>(cnt_s_vec[s]);
        if (lhs > cnt) {
            h_mu = 0;
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

    for (auto& v: candidates.at(model::ComplexLocation::PATCH_SURF)) {
        cmult *= v.second.rateMult(pdef->complexStates(v.first));
    }

    // Volume part, if applicable
    if (defcsr().inside() or defcsr().outside()) {
        const auto& lhs_vec = defcsr().inside() ? pdef->complexsreac_lhs_I(lidx)
                                                : pdef->complexsreac_lhs_O(lidx);
        const auto& cnt_vec = defcsr().inside() ? pPatch.iComp()->def()->pools()
                                                : pPatch.oComp()->def()->pools();
        if (h_mu > 0 and cmult > 0) {
            for (auto s: lhs_vec.range()) {
                const uint lhs = lhs_vec[s];
                if (lhs == 0) {
                    continue;
                }
                const uint cnt = cnt_vec[s];
                if (lhs > cnt) {
                    h_mu = 0;
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
                }
                }
            }
        }

        const auto loc = defcsr().inside() ? model::ComplexLocation::PATCH_IN
                                           : model::ComplexLocation::PATCH_OUT;
        const auto cdef = defcsr().inside() ? pPatch.iComp()->def() : pPatch.oComp()->def();
        for (auto& v: candidates.at(loc)) {
            cmult *= v.second.rateMult(cdef->complexStates(v.first));
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * cmult * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& ComplexSReac::apply() {
    solver::Patchdef* pdef = pPatch.def();
    Comp* icomp = pPatch.iComp();
    Comp* ocomp = pPatch.oComp();
    solver::complexsreac_local_id lidx = pdef->complexsreacG2L(defcsr().gidx());
    auto& rng = pdef->statedef().rng();

    // Surface species
    const auto& upd_s_vec = pdef->complexsreac_upd_S(lidx);
    const auto& cnt_s_vec = pdef->pools();
    for (auto s: upd_s_vec.range()) {
        if (pdef->clamped(s)) {
            continue;
        }
        const int upd = upd_s_vec[s];
        if (upd == 0) {
            continue;
        }
        const int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        AssertLog(nc >= 0);
        pdef->setCount(s, static_cast<double>(nc));
    }

    // Surface complexes
    for (auto& v: candidates.at(model::ComplexLocation::PATCH_SURF)) {
        for (auto& event: v.second.selectEvents(rng)) {
            if (event.first->type() == solver::UPDEvent) {
                const auto ev = std::dynamic_pointer_cast<const solver::ComplexUpdateEventdef>(
                    event.first);
                // Updates
                const auto& state = pdef->complexStates(v.first).at(event.second);
                pdef->updateComplex(v.first, event.second, ev->getUpdate(state, rng));

                switch (ev->destLoc()) {
                case model::ComplexLocation::PATCH_IN: {
                    AssertLog(icomp != nullptr);
                    icomp->def()->addComplex(v.first,
                                             event.second,
                                             pdef->complexStates(v.first).at(event.second));
                    pdef->removeComplex(v.first, event.second);
                    break;
                }
                case model::ComplexLocation::PATCH_OUT: {
                    AssertLog(ocomp != nullptr);
                    ocomp->def()->addComplex(v.first,
                                             event.second,
                                             pdef->complexStates(v.first).at(event.second));
                    pdef->removeComplex(v.first, event.second);
                    break;
                }
                default:
                    break;
                }
            } else {
                // Deletions
                pdef->removeComplex(v.first, event.second);
            }
        }
    }
    // Creations
    for (auto& ce: defcsr().creEvents(model::ComplexLocation::PATCH_SURF)) {
        pdef->addComplex(ce->complexIdx(), pPatch.solver()->getNextComplexInd(), ce->init());
    }

    // Volume part
    if (icomp != nullptr) {
        const auto& upd_vec = pdef->complexsreac_upd_I(lidx);
        const auto& cnt_vec = icomp->def()->pools();
        for (auto s: upd_vec.range()) {
            if (icomp->def()->clamped(s)) {
                continue;
            }
            const int upd = upd_vec[s];
            if (upd == 0) {
                continue;
            }
            const int nc = static_cast<int>(cnt_vec[s]) + upd;
            AssertLog(nc >= 0);
            icomp->def()->setCount(s, static_cast<double>(nc));
        }
        for (auto& ce: defcsr().creEvents(model::ComplexLocation::PATCH_IN)) {
            icomp->def()->addComplex(ce->complexIdx(),
                                     pPatch.solver()->getNextComplexInd(),
                                     ce->init());
        }
    }
    if (ocomp != nullptr) {
        const auto& upd_vec = pdef->complexsreac_upd_O(lidx);
        const auto& cnt_vec = ocomp->def()->pools();
        for (auto s: upd_vec.range()) {
            if (ocomp->def()->clamped(s)) {
                continue;
            }
            const int upd = upd_vec[s];
            if (upd == 0) {
                continue;
            }
            const int nc = static_cast<int>(cnt_vec[s]) + upd;
            AssertLog(nc >= 0);
            ocomp->def()->setCount(s, static_cast<double>(nc));
        }
        for (auto& ce: defcsr().creEvents(model::ComplexLocation::PATCH_OUT)) {
            ocomp->def()->addComplex(ce->complexIdx(),
                                     pPatch.solver()->getNextComplexInd(),
                                     ce->init());
        }
    }
    // Volume LHS complexes
    if (defcsr().inside() or defcsr().outside()) {
        const auto comp = defcsr().inside() ? pPatch.iComp() : pPatch.oComp();

        const auto loc = defcsr().inside() ? model::ComplexLocation::PATCH_IN
                                           : model::ComplexLocation::PATCH_OUT;
        for (auto& v: candidates.at(loc)) {
            for (auto& event: v.second.selectEvents(rng)) {
                if (event.first->type() == solver::UPDEvent) {
                    const auto ev = std::dynamic_pointer_cast<const solver::ComplexUpdateEventdef>(
                        event.first);
                    // Updates
                    const auto& state = comp->def()->complexStates(v.first).at(event.second);
                    comp->def()->updateComplex(v.first, event.second, ev->getUpdate(state, rng));

                    if (ev->destLoc() != loc) {
                        switch (ev->destLoc()) {
                        case model::ComplexLocation::PATCH_SURF: {
                            pdef->addComplex(v.first,
                                             event.second,
                                             comp->def()->complexStates(v.first).at(event.second));
                            comp->def()->removeComplex(v.first, event.second);
                            break;
                        }
                        case model::ComplexLocation::PATCH_IN:
                        case model::ComplexLocation::PATCH_OUT: {
                            const auto other = defcsr().inside() ? pPatch.oComp() : pPatch.iComp();
                            AssertLog(other != nullptr);
                            other->def()->addComplex(v.first,
                                                     event.second,
                                                     comp->def()->complexStates(v.first).at(
                                                         event.second));
                            comp->def()->removeComplex(v.first, event.second);
                            break;
                        }
                        default:
                            break;
                        }
                    }
                } else {
                    // Deletions
                    comp->def()->removeComplex(v.first, event.second);
                }
            }
        }
    }

    rExtent++;
    return pUpdVec;
}

}  // namespace steps::wmdirect
