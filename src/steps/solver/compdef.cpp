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

#include "compdef.hpp"

#include "compdef.hpp"
#include "complexdef.hpp"
#include "complexreacdef.hpp"
#include "diffdef.hpp"
#include "geom/comp.hpp"
#include "model/diff.hpp"
#include "model/model.hpp"
#include "model/volsys.hpp"
#include "patchdef.hpp"
#include "reacdef.hpp"
#include "statedef.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

namespace steps::solver {
Compdef::Compdef(Statedef& sd, comp_global_id idx, wm::Comp& c)
    : pStatedef(sd)
    , pName(c.getID())
    , pVol(c.getVol())
    , pIdx(idx)
    , pCvsys(c.getVolsys()) {
    uint nspecs = pStatedef.countSpecs();
    pSpec_G2L.container().resize(nspecs);

    uint nreacs = pStatedef.countReacs();
    pReac_G2L.container().resize(nreacs);

    uint ncreacs = pStatedef.countComplexReacs();
    pComplexReac_G2L.container().resize(ncreacs);

    uint ndiffs = pStatedef.countDiffs();
    pDiff_G2L.container().resize(ndiffs);

    uint nvesbinds = pStatedef.countVesBinds();
    pVesBind_G2L.container().resize(nvesbinds);

    uint nvesunbinds = pStatedef.countVesUnbinds();
    pVesUnbind_G2L.container().resize(nvesunbinds);

    // We can already setup the complexes because they do not use local indices
    uint nbComplexes = pStatedef.countComplexes();
    pComplexStates.container().resize(nbComplexes);
    pFilters.container().resize(nbComplexes);
    pFiltersMap.container().resize(nbComplexes);
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::checkpoint(std::fstream& cp_file) const {
    util::checkpoint(cp_file, pVol);
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pReacKcst);
    util::checkpoint(cp_file, pReacFlags);
    util::checkpoint(cp_file, pComplexReacKcst);
    util::checkpoint(cp_file, pComplexReacFlags);
    util::checkpoint(cp_file, pDiffDcst);
    util::checkpoint(cp_file, pComplexStates);
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::restore(std::fstream& cp_file) {
    util::restore(cp_file, pVol);
    util::restore(cp_file, pPoolCount);
    util::restore(cp_file, pPoolFlags);
    util::restore(cp_file, pReacKcst);
    util::restore(cp_file, pReacFlags);
    util::restore(cp_file, pComplexReacKcst);
    util::restore(cp_file, pComplexReacFlags);
    util::restore(cp_file, pDiffDcst);
    util::restore(cp_file, pComplexStates);
    for (auto cid: pFilters.range()) {
        pFilters[cid].container().clear();
        pFiltersMap[cid].clear();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setup_references() {
    AssertLog(pSetupRefsdone == false);
    AssertLog(pSetupIndsdone == false);

    const uint ngspecs = pStatedef.countSpecs();
    const uint ngreacs = pStatedef.countReacs();
    const uint ngcreacs = pStatedef.countComplexReacs();
    const uint ngdiffs = pStatedef.countDiffs();
    const uint ngvesbinds = pStatedef.countVesBinds();
    const uint ngvesunbinds = pStatedef.countVesUnbinds();

    if (ngspecs == 0) {
        AssertLog(pSpec_G2L.container().empty());
    }
    if (ngreacs == 0) {
        AssertLog(pReac_G2L.container().empty());
    }
    if (ngdiffs == 0) {
        AssertLog(pDiff_G2L.container().empty());
    }
    if (ngvesbinds == 0) {
        AssertLog(pVesBind_G2L.container().empty());
    }
    if (ngvesunbinds == 0) {
        AssertLog(pVesUnbind_G2L.container().empty());
    }

    // Importantly  assumes that all species from patch sreacs have
    // been added first. Statedef calls setup on patches, which add Specs to their
    // inner and outer compartments.

    // set up local reac indices ////vsys also has _countReacs and Reac *
    // _getReac(lidx) The local reacs index count (pReacsN) starts at 0.
    for (auto const& v: pCvsys) {
        const auto& vreacs = pStatedef.model().getVolsys(v)._getAllReacs();
        if (ngreacs == 0) {
            AssertLog(vreacs.empty() == true);
        }
        for (auto const& [_, reac]: vreacs) {
            reac_global_id gidx = pStatedef.getReacIdx(*reac);
            AssertLog(gidx < ngreacs);
            if (reacG2L(gidx).valid()) {
                continue;
            }
            pReac_G2L[gidx] = reac_local_id(pReacsN++);
        }
        const auto& vdiffs = pStatedef.model().getVolsys(v)._getAllDiffs();
        if (ngdiffs == 0) {
            AssertLog(vdiffs.empty() == true);
        }
        for (auto const& [_, diff]: vdiffs) {
            diff_global_id gidx = pStatedef.getDiffIdx(*diff);
            AssertLog(gidx < ngdiffs);
            if (diffG2L(gidx).valid()) {
                continue;
            }
            pDiff_G2L[gidx] = diff_local_id(pDiffsN++);
        }

        const auto& vvesbinds = pStatedef.model().getVolsys(v)._getAllVesBinds();
        if (ngvesbinds == 0) {
            AssertLog(vvesbinds.empty() == true);
        }
        for (auto const& [_, vb]: vvesbinds) {
            vesbind_global_id gidx = pStatedef.getVesBindIdx(*vb);
            AssertLog(gidx < ngvesbinds);
            if (vesbindG2L(gidx).valid()) {
                continue;
            }
            pVesBind_G2L[gidx] = vesbind_local_id(pVesBindsN++);
        }

        const auto& vvesunbinds = pStatedef.model().getVolsys(v)._getAllVesUnbinds();
        if (ngvesunbinds == 0) {
            AssertLog(vvesunbinds.empty() == true);
        }
        for (const auto& [_, vub]: vvesunbinds) {
            vesunbind_global_id gidx = pStatedef.getVesUnbindIdx(*vub);
            AssertLog(gidx < ngvesunbinds);
            if (vesunbindG2L(gidx).valid()) {
                continue;
            }
            pVesUnbind_G2L[gidx] = vesunbind_local_id(pVesUnbindsN++);
        }

        auto vcreacs = pStatedef.model().getVolsys(v)._getAllComplexReacs();
        if (ngcreacs == 0) {
            AssertLog(vcreacs.empty());
        }
        for (auto const& r: vcreacs) {
            complexreac_global_id gidx = pStatedef.getComplexReacIdx(r.second);
            AssertLog(gidx < ngcreacs);
            if (pComplexReac_G2L[gidx].valid()) {
                continue;
            }
            pComplexReac_G2L[gidx] = complexreac_local_id(pComplexReacsN++);
        }
    }

    // now add all species that appear in all reactions, diffusions that
    // can occur in this compartment
    // NOTE: Patchdef setups have called addSpec() to already add some
    // species (from surface reactions)
    for (auto r: reac_global_id::range(ngreacs)) {
        if (reacG2L(r).unknown()) {
            continue;
        }

        const Reacdef& rdef = pStatedef.reacdef(r);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (rdef.reqspec(s) == true) {
                addSpec(s);
            }
        }
    }
    for (auto d: diff_global_id::range(ngdiffs)) {
        if (diffG2L(d).unknown()) {
            continue;
        }

        const Diffdef& ddef = pStatedef.diffdef(d);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (ddef.reqspec(s)) {
                addSpec(s);
            }
        }
    }
    for (auto r: complexreac_global_id::range(ngcreacs)) {
        if (pComplexReac_G2L[r].unknown()) {
            continue;
        }
        const ComplexReacdef& crdef = pStatedef.complexreacdef(r);
        for (auto s: spec_global_id::range(ngspecs)) {
            if (crdef.reqspec(s)) {
                addSpec(s);
            }
        }
    }

    // Complexes were already added because we do not use local indices

    pSetupRefsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setup_indices() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == false);

    uint ngspecs = pStatedef.countSpecs();
    uint ngreacs = pStatedef.countReacs();
    uint ngcreacs = pStatedef.countComplexReacs();
    uint ngdiffs = pStatedef.countDiffs();
    uint ngvesbinds = pStatedef.countVesBinds();
    uint ngvesunbinds = pStatedef.countVesUnbinds();

    // Set up local indices
    if (countSpecs() != 0) {
        pSpec_L2G.container().resize(countSpecs());
        for (auto i: spec_global_id::range(ngspecs)) {
            spec_local_id lidx = specG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pSpec_L2G[lidx] = i;
        }
    }

    if (countReacs() != 0) {
        pReac_L2G.container().resize(countReacs());
        for (auto i: reac_global_id::range(ngreacs)) {
            reac_local_id lidx = reacG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pReac_L2G[lidx] = i;
        }
        uint arrsize = countSpecs() * countReacs();
        pReac_DEP_Spec.resize(arrsize);
        pReac_LHS_Spec.resize(arrsize);
        pReac_UPD_Spec.resize(arrsize);
        for (auto ri: reac_local_id::range(countReacs())) {
            const Reacdef& rdef = reacdef(ri);
            for (auto si: spec_global_id::range(ngspecs)) {
                if (rdef.reqspec(si) == false) {
                    continue;
                }
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());

                uint aridx = _IDX_Reac_Spec(ri, sil);
                pReac_DEP_Spec[aridx] = rdef.dep(si);
                pReac_LHS_Spec[aridx] = rdef.lhs(si);
                pReac_UPD_Spec[aridx] = rdef.upd(si);
            }
        }
    }

    if (countComplexReacs() != 0) {
        pComplexReac_L2G.container().resize(countComplexReacs());
        for (auto i: complexreac_global_id::range(ngcreacs)) {
            complexreac_local_id lidx = pComplexReac_G2L[i];
            if (lidx.unknown()) {
                continue;
            }
            pComplexReac_L2G[lidx] = i;
        }
        uint arrsize = countSpecs() * pComplexReacsN;
        pComplexReac_DEP_Spec.resize(arrsize);
        pComplexReac_LHS_Spec.resize(arrsize);
        pComplexReac_UPD_Spec.resize(arrsize);
        for (auto cri: complexreac_local_id::range(countComplexReacs())) {
            ComplexReacdef& rdef = complexreacdef(cri);
            for (auto si: spec_global_id::range(ngspecs)) {
                if (rdef.reqspec(si) == false) {
                    continue;
                }
                spec_local_id sil = pSpec_G2L[si];
                AssertLog(sil.valid());

                uint aridx = _IDX_ComplexReac_Spec(cri, sil);
                pComplexReac_DEP_Spec[aridx] = rdef.dep(si);
                pComplexReac_LHS_Spec[aridx] = rdef.lhs(si);
                pComplexReac_UPD_Spec[aridx] = rdef.upd(si);
            }
        }
    }

    if (countDiffs() != 0) {
        pDiff_L2G.container().resize(countDiffs());
        for (auto i: diff_global_id::range(ngdiffs)) {
            diff_local_id lidx = diffG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pDiff_L2G[lidx] = i;
        }

        uint arrsize = countSpecs() * countDiffs();
        pDiff_DEP_Spec.resize(arrsize);
        pDiff_LIG.container().resize(countDiffs());
        for (auto di: diff_local_id::range(countDiffs())) {
            const Diffdef& ddef = diffdef(di);
            pDiff_LIG[di] = specG2L(ddef.lig());
            for (auto si: spec_global_id::range(ngspecs)) {
                if (ddef.reqspec(si) == false) {
                    continue;
                }
                spec_local_id sil = specG2L(si);
                AssertLog(sil.valid());
                uint aridx = _IDX_Diff_Spec(di, sil);
                pDiff_DEP_Spec[aridx] = ddef.dep(si);
            }
        }
    }

    if (countVesBinds() != 0) {
        pVesBind_L2G.container().resize(countVesBinds());
        for (auto i: vesbind_global_id::range(ngvesbinds)) {
            vesbind_local_id lidx = vesbindG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pVesBind_L2G[lidx] = i;
        }
    }

    if (countVesUnbinds() != 0) {
        pVesUnbind_L2G.container().resize(countVesUnbinds());
        for (auto i: vesunbind_global_id::range(ngvesunbinds)) {
            vesunbind_local_id lidx = vesunbindG2L(i);
            if (lidx.unknown()) {
                continue;
            }
            pVesUnbind_L2G[lidx] = i;
        }
    }

    if (countComplexReacs() != 0) {
        pComplexReacFlags.container().resize(pComplexReacsN);

        // Finally initialise constants to user-supplied values
        pComplexReacKcst.container().resize(pComplexReacsN);

        for (auto i: complexreac_local_id::range(pComplexReacsN)) {
            ComplexReacdef& reac = complexreacdef(i);
            pComplexReacKcst[i] = reac.kcst();
        }
    }

    // Initialise the pools and flags members to zeros.
    if (countSpecs() != 0) {
        pPoolCount.container().resize(countSpecs());
        pPoolFlags.container().resize(countSpecs());
    }
    if (countReacs() != 0) {
        pReacFlags.container().resize(countReacs());

        // Finally initialise constants to user-supplied values
        pReacKcst.container().resize(countReacs());

        for (auto i: reac_local_id::range(countReacs())) {
            // reacdef() returns global Reacdef by local index
            pReacKcst[i] = reacdef(i).kcst();
        }
    }

    if (countDiffs() != 0) {
        pDiffDcst.container().resize(countDiffs());

        for (auto i: diff_local_id::range(countDiffs())) {
            // diffdef() returns global Diffdef by local index
            pDiffDcst[i] = diffdef(i).dcst();
        }
    }
    pSetupIndsdone = true;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::addIPatchdef(Patchdef& p) {
    // Make some checks.
    AssertLog(p.ocompdef() == this);
    // Check whether it's already included.
    auto ip_end = pIPatches.end();
    if (std::find(pIPatches.begin(), ip_end, &p) != ip_end) {
        return;
    }
    auto op_end = pOPatches.end();
    AssertLog(std::find(pOPatches.begin(), op_end, &p) == op_end);
    // Include.
    pIPatches.push_back(&p);
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::addOPatchdef(Patchdef& p) {
    // Make some checks.
    AssertLog(p.icompdef() == this);
    // Check whether it's already included.
    auto op_end = pOPatches.end();
    if (std::find(pOPatches.begin(), op_end, &p) != op_end) {
        return;
    }
    auto ip_end = pIPatches.end();
    AssertLog(std::find(pIPatches.begin(), ip_end, &p) == ip_end);
    // Include.
    pOPatches.push_back(&p);
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::addSpec(spec_global_id gidx) {
    AssertLog(pSetupIndsdone == false);
    if (specG2L(gidx).valid()) {
        return;
    }
    pSpec_G2L[gidx] = spec_local_id(pSpecsN++);
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setVol(double v) {
    AssertLog(v > 0.0);
    pVol = v;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::reset() {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0.0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);
    std::fill(pReacFlags.begin(), pReacFlags.end(), 0);
    std::fill(pComplexReacFlags.begin(), pComplexReacFlags.end(), 0);
    for (auto& map: pComplexStates) {
        map.clear();
    }
    for (auto& filts: pFilters) {
        for (auto filt: filts) {
            filt->reset();
        }
    }
    for (auto i: reac_local_id::range(countReacs())) {
        pReacKcst[i] = reacdef(i).kcst();
    }
    for (auto i: complexreac_local_id::range(pComplexReacsN)) {
        ComplexReacdef& creac = complexreacdef(i);
        pComplexReacKcst[i] = creac.kcst();
    }
    for (auto i: diff_local_id::range(countDiffs())) {
        pDiffDcst[i] = diffdef(i).dcst();
    }
}

////////////////////////////////////////////////////////////////////////////////

std::shared_ptr<ComplexFilter> Compdef::GetFilter(
    const complex_global_id& cmplIdx,
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>&
        filts) {
    auto it = pFiltersMap[cmplIdx].find(filts);
    if (it == pFiltersMap[cmplIdx].end()) {
        complex_filter_id fid(pFilters[cmplIdx].size());
        pFilters[cmplIdx].container().emplace_back(
            new ComplexFilter(filts, fid, pComplexStates[cmplIdx]));
        it = pFiltersMap[cmplIdx].emplace(filts, pFilters[cmplIdx].back()).first;
    }
    return it->second;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::updateComplex(complex_global_id cmplIdx,
                            complex_individual_id stIdx,
                            const util::strongid_vector<complex_substate_id, int>& upd) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(cmplIdx.get() < pComplexStates.container().size());
    AssertLog(pComplexStates[cmplIdx].count(stIdx) > 0);

    auto& state = pComplexStates[cmplIdx].at(stIdx);
    std::vector<complex_substate_id> modifInds;
    for (auto sus: upd.range()) {
        if (upd[sus] != 0) {
            modifInds.push_back(sus);
            state[sus] += upd[sus];
        }
    }

    if (not modifInds.empty()) {
        for (auto filt: pFilters[cmplIdx]) {
            for (auto& sus: modifInds) {
                if (filt->dependsOnSus(sus) or filt->matchAll()) {
                    filt->toUpdate(state.ind());
                    break;
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::removeComplex(complex_global_id cmplIdx, complex_individual_id stIdx) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(cmplIdx.get() < pComplexStates.container().size());
    AssertLog(pComplexStates[cmplIdx].count(stIdx) > 0);

    // Empty the vector at this position and add the position to the list of holes
    pComplexStates[cmplIdx].erase(stIdx);

    for (auto filt: pFilters[cmplIdx]) {
        filt->toUpdate(stIdx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::addComplex(complex_global_id cmplIdx,
                         complex_individual_id stIdx,
                         const util::strongid_vector<complex_substate_id, uint>& init) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(cmplIdx.get() < pComplexStates.container().size());

    pComplexStates[cmplIdx].emplace(stIdx, ComplexState(init, stIdx));

    for (auto filt: pFilters[cmplIdx]) {
        filt->toUpdate(stIdx);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::clearComplexFilterUpdates() {
    for (auto& filters: pFilters) {
        for (auto filt: filters) {
            filt->clearLastUpdates();
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setCount(spec_local_id slidx, double count) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(slidx < countSpecs());
    AssertLog(count >= 0.0);
    pPoolCount[slidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setClamped(spec_local_id slidx, bool clamp) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(slidx < countSpecs());
    if (clamp == true) {
        pPoolFlags[slidx] |= CLAMPED;
    } else {
        pPoolFlags[slidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

strongid_span<spec_local_id, const uint> Compdef::reac_lhs(reac_local_id rlidx) const {
    AssertLog(rlidx < countReacs());
    return to_strong_span<spec_local_id>(pReac_LHS_Spec, rlidx, countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

strongid_span<spec_local_id, const int> Compdef::reac_upd(reac_local_id rlidx) const {
    AssertLog(rlidx < countReacs());
    return to_strong_span<spec_local_id>(pReac_UPD_Spec, rlidx, countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int Compdef::reac_dep(reac_local_id rlidx, spec_local_id slidx) const {
    return pReac_DEP_Spec[slidx.get() + (rlidx.get() * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

uint Compdef::diff_dep(diff_local_id dlidx, spec_local_id slidx) const {
    return pDiff_DEP_Spec[slidx.get() + (dlidx.get() * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

Reacdef& Compdef::reacdef(reac_local_id rlidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(rlidx < countReacs());
    return pStatedef.reacdef(reacL2G(rlidx));
}

////////////////////////////////////////////////////////////////////////////////

ComplexReacdef& Compdef::complexreacdef(complexreac_local_id rlidx) const {
    AssertLog(pSetupRefsdone);
    AssertLog(rlidx < countComplexReacs());
    return pStatedef.complexreacdef(pComplexReac_L2G[rlidx]);
}

////////////////////////////////////////////////////////////////////////////////

strongid_span<spec_local_id, const uint> Compdef::complexreac_lhs(
    complexreac_local_id rlidx) const {
    AssertLog(rlidx < countComplexReacs());
    return to_strong_span<spec_local_id>(pComplexReac_LHS_Spec, rlidx, countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

strongid_span<spec_local_id, const int> Compdef::complexreac_upd(complexreac_local_id rlidx) const {
    AssertLog(rlidx < countComplexReacs());
    return to_strong_span<spec_local_id>(pComplexReac_UPD_Spec, rlidx, countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int Compdef::complexreac_dep(complexreac_local_id rlidx, spec_local_id slidx) const {
    return pComplexReac_DEP_Spec[slidx.get() + (rlidx.get() * pSpecsN)];
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setKcst(reac_local_id rlidx, double kcst) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(rlidx < countReacs());
    AssertLog(kcst >= 0.0);
    pReacKcst[rlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setDcst(diff_local_id dlidx, double dcst) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(dlidx < countDiffs());
    AssertLog(dcst >= 0.0);
    pDiffDcst[dlidx] = dcst;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setActive(reac_local_id rlidx, bool active) {
    AssertLog(pSetupRefsdone == true);
    AssertLog(pSetupIndsdone == true);
    AssertLog(rlidx < countReacs());
    if (active == true) {
        pReacFlags[rlidx] &= ~INACTIVATED;
    } else {
        pReacFlags[rlidx] |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setComplexReacKcst(complexreac_local_id rlidx, double kcst) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(rlidx < pComplexReacsN);
    AssertLog(kcst >= 0.0);
    pComplexReacKcst[rlidx] = kcst;
}

////////////////////////////////////////////////////////////////////////////////

void Compdef::setComplexReacActive(complexreac_local_id rlidx, bool active) {
    AssertLog(pSetupRefsdone);
    AssertLog(pSetupIndsdone);
    AssertLog(rlidx < countComplexReacs());
    if (active) {
        pComplexReacFlags[rlidx] &= ~INACTIVATED;
    } else {
        pComplexReacFlags[rlidx] |= INACTIVATED;
    }
}

////////////////////////////////////////////////////////////////////////////////

Diffdef& Compdef::diffdef(diff_local_id dlidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(dlidx < countDiffs());
    return pStatedef.diffdef(diffL2G(dlidx));
}

////////////////////////////////////////////////////////////////////////////////

VesBinddef& Compdef::vesbinddef(vesbind_local_id vblidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(vblidx < countVesBinds());
    return pStatedef.vesbinddef(vesbindL2G(vblidx));
}

////////////////////////////////////////////////////////////////////////////////

VesUnbinddef& Compdef::vesunbinddef(vesunbind_local_id vblidx) const {
    AssertLog(pSetupRefsdone == true);
    AssertLog(vblidx < countVesUnbinds());
    return pStatedef.vesunbinddef(vesunbindL2G(vblidx));
}

}  // namespace steps::solver
