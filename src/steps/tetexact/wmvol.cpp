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
#include "wmvol.hpp"

#include "diff.hpp"
#include "math/constants.hpp"
#include "reac.hpp"
#include "solver/reacdef.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "tri.hpp"
#include "util/checkpointing.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetexact {


////////////////////////////////////////////////////////////////////////////////

WmVol::WmVol(tetrahedron_global_id idx, solver::Compdef* cdef, double vol)
    : pIdx(idx)
    , pCompdef(cdef)
    , pVol(vol) {
    AssertLog(pCompdef != nullptr);
    AssertLog(pVol > 0.0);

    // Based on compartment definition, build other structures.
    auto nspecs = compdef()->countSpecs();
    pPoolCount.container().resize(nspecs, 0);
    pPoolFlags.container().resize(nspecs, 0);
    pKProcs.resize(compdef()->countReacs());
}

////////////////////////////////////////////////////////////////////////////////

WmVol::~WmVol() {
    // Delete reaction rules.
    for (auto const& i: pKProcs) {
        delete i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pPoolFlags);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolCount);
    util::restore(cp_file, pPoolFlags);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setNextTri(Tri* t) {
    pNextTris.push_back(t);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setupKProcs(Tetexact* tex) {
    uint j = 0;

    // Note: ignoring diffusion KProcs

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (auto i: solver::reac_local_id::range(nreacs)) {
        auto& rdef = compdef()->reacdef(i);
        auto* r = new Reac(&rdef, this);
        pKProcs[j++] = r;
        tex->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::reset() {
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);

    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }
}

////////////////////////////////////////////////////////////////////////////////

double WmVol::conc(solver::spec_global_id gidx) const {
    solver::spec_local_id lspidx = compdef()->specG2L(gidx);
    double n = pPoolCount[lspidx];
    return n / (1.0e3 * pVol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setCount(solver::spec_local_id lidx, uint count) {
    pPoolCount.at(lidx) = count;
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::incCount(solver::spec_local_id lidx, int inc) {
    AssertLog(lidx < compdef()->countSpecs());
#ifndef NDEBUG
    uint old_count = pPoolCount[lidx];
#endif

    pPoolCount[lidx] += inc;

#ifndef NDEBUG
    uint new_count = pPoolCount[lidx];
    AssertLog((inc >= 0 && new_count >= old_count) || (inc < 0 && new_count < old_count));
#endif
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setClamped(solver::spec_local_id lidx, bool clamp) {
    if (clamp) {
        pPoolFlags[lidx] |= CLAMPED;
    } else {
        pPoolFlags[lidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

Reac& WmVol::reac(solver::reac_local_id lidx) const {
    return *dynamic_cast<Reac*>(pKProcs.at(lidx.get()));
}

}  // namespace steps::tetexact
