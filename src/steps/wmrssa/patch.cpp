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
#include "patch.hpp"
#include "solver/statedef.hpp"
#include "wmrssa.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::wmrssa {

////////////////////////////////////////////////////////////////////////////////

Patch::Patch(solver::Patchdef* patchdef, Comp* icomp, Comp* ocomp)
    : pPatchdef(patchdef)
    , pIComp(icomp)
    , pOComp(ocomp) {
    AssertLog(pPatchdef != nullptr);
    if (iComp() != nullptr) {
        iComp()->addIPatch(this);
    }
    if (oComp() != nullptr) {
        oComp()->addOPatch(this);
    }
    uint nspecs = patchdef->countSpecs();
    pPoolLB.container().resize(nspecs);
    pPoolUB.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch() {
    for (auto const& k: pKProcs) {
        delete k;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patch::checkpoint(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->checkpoint(cp_file);
    }
    util::checkpoint(cp_file, pPoolLB);
    util::checkpoint(cp_file, pPoolUB);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::restore(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->restore(cp_file);
    }
    util::restore(cp_file, pPoolLB);
    util::restore(cp_file, pPoolUB);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setupKProcs(Wmrssa* wmd) {
    // Create surface reaction kproc's.
    uint nsreacs = def()->countSReacs();
    pKProcs.resize(nsreacs);
    for (auto i: solver::sreac_local_id::range(nsreacs)) {
        auto& srdef = def()->sreacdef(i);
        auto* sr = new SReac(&srdef, this);
        pKProcs[i.get()] = sr;
        wmd->addKProc(sr);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setupDeps() {
    for (auto const& kproc: pKProcs) {
        kproc->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////

KProc* Patch::sreac(solver::sreac_local_id lsridx) const {
    assert(lsridx.get() < pKProcs.size());
    return pKProcs[lsridx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Patch::reset() {
    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setBounds(solver::spec_local_id i, int nc) {
    constexpr double delta = .05;
    if (nc > 3 / delta) {
        pPoolLB[i] = nc * (1 - delta);
        pPoolUB[i] = nc * (1 + delta);
    } else if (nc > 3) {
        pPoolLB[i] = nc - 3;
        pPoolUB[i] = nc + 3;
    } else if (nc > 0) {
        pPoolLB[i] = 1;
        pPoolUB[i] = 2 * nc;
    } else {
        pPoolLB[i] = 0;
        pPoolUB[i] = 0;
    }
    pPoolLB[i] -= delta;
    pPoolUB[i] += delta;
    /*pPoolLB[i] = std::max(std::min(nc*(1 - delta), nc - 3.), 0.); //nc/(3 -
    delta - 2*(1 - delta)/(1 + 2./(nc+1))); // pPoolUB[i] = std::max(nc*(1 +
    delta), nc + 3.); //nc*(3 - delta - 2*(1 - delta)/(1 + 2./(nc+1))); /*/
}

////////////////////////////////////////////////////////////////////////////////

bool Patch::isOutOfBound(solver::spec_local_id i, int nc) {
    AssertLog(i < def()->countSpecs());
    if (nc > pPoolLB[i] && nc < pPoolUB[i]) {
        return false;
    }
    setBounds(i, nc);
    return true;
}

////////////////////////////////////////////////////////////////////////////////

const util::strongid_vector<solver::spec_local_id, double>& Patch::pools(
    wmrssa::PropensityRSSA prssa) const {
    switch (prssa) {
    case wmrssa::CURRENT:
        return def()->pools();
    case wmrssa::LOWERBOUND:
        return pPoolLB;
    case wmrssa::BOUNDS:
        return pPoolUB;
    default:
        AssertLog(false);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setupSpecDeps() {
    uint nspecs = def()->countSpecs();
    localSpecUpdKProcs.resize(nspecs);
    for (auto slidx: solver::spec_local_id::range(nspecs)) {
        solver::spec_global_id sgidx = def()->specL2G(slidx);
        for (auto const& k: pKProcs) {
            if (k->depSpecPatch(sgidx, this)) {
                localSpecUpdKProcs[slidx.get()].push_back(k);
            }
        }
        if (pIComp != nullptr) {
            for (auto const& k: pIComp->kprocs()) {
                if (k->depSpecPatch(sgidx, this)) {
                    localSpecUpdKProcs[slidx.get()].push_back(k);
                }
            }
        }
        if (pOComp != nullptr) {
            for (auto const& k: pOComp->kprocs()) {
                if (k->depSpecPatch(sgidx, this)) {
                    localSpecUpdKProcs[slidx.get()].push_back(k);
                }
            }
        }
    }
}

}  // namespace steps::wmrssa
