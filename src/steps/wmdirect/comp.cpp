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
#include "comp.hpp"
#include "complexreac.hpp"
#include "reac.hpp"
#include "wmdirect.hpp"
// logging
#include "util/error.hpp"

namespace steps::wmdirect {

Comp::Comp(solver::Compdef* compdef, Wmdirect* solver)
    : pSolver(solver)
    , pCompdef(compdef) {
    AssertLog(pCompdef != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp() {
    for (auto const& k: pKProcs) {
        delete k;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::checkpoint(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->checkpoint(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::restore(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::reset() {
    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setupKProcs(Wmdirect* wmd) {
    // Create reaction kproc's.
    uint nreacs = def()->countReacs();
    uint ncreacs = def()->countComplexReacs();
    pKProcs.resize(nreacs + ncreacs);
    for (auto i: solver::reac_local_id::range(nreacs)) {
        solver::Reacdef& rdef = def()->reacdef(i);
        auto* r = new Reac(&rdef, this);
        pKProcs[i.get()] = r;
        wmd->addKProc(r);
    }

    for (auto i: solver::complexreac_local_id::range(ncreacs)) {
        solver::ComplexReacdef& rdef = def()->complexreacdef(i);
        auto* r = new ComplexReac(rdef, *this);
        pKProcs[nreacs + i.get()] = r;
        wmd->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

KProc* Comp::reac(solver::reac_local_id lridx) const {
    AssertLog(lridx.get() < pKProcs.size());
    return pKProcs[lridx.get()];
}

////////////////////////////////////////////////////////////////////////////////

KProc* Comp::reac(solver::complexreac_local_id lridx) const {
    uint idx = def()->countReacs() + lridx.get();
    AssertLog(idx < pKProcs.size());
    return pKProcs[idx];
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setupDeps() {
    for (auto const& kproc: pKProcs) {
        kproc->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addIPatch(Patch* p) {
    AssertLog(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addOPatch(Patch* p) {
    AssertLog(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}

}  // namespace steps::wmdirect
