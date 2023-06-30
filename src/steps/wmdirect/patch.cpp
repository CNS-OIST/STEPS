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
#include "wmdirect.hpp"
// logging
#include "util/error.hpp"
#include <easylogging++.h>

namespace steps::wmdirect {

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
}

////////////////////////////////////////////////////////////////////////////////

void Patch::restore(std::fstream& cp_file) {
    for (auto const& k: pKProcs) {
        k->restore(cp_file);
    }
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setupKProcs(Wmdirect* wmd) {
    // Create surface reaction kproc's.
    uint nsreacs = def()->countSReacs();
    pKProcs.resize(nsreacs);
    for (auto i: solver::sreac_local_id::range(nsreacs)) {
        solver::SReacdef* srdef = def()->sreacdef(i);
        auto* sr = new SReac(srdef, this);
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
    AssertLog(lsridx.get() < pKProcs.size());
    return pKProcs[lsridx.get()];
}

////////////////////////////////////////////////////////////////////////////////

void Patch::reset() {
    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }
}

}  // namespace steps::wmdirect
