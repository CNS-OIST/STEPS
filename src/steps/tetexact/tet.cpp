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

#include "tet.hpp"

#include "diff.hpp"
#include "reac.hpp"
#include "solver/diffdef.hpp"
#include "tetexact.hpp"
#include "tri.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetexact {

Tet::Tet(tetrahedron_global_id idx,
         solver::Compdef* cdef,
         double vol,
         double a0,
         double a1,
         double a2,
         double a3,
         double d0,
         double d1,
         double d2,
         double d3,
         tetrahedron_global_id tet0,
         tetrahedron_global_id tet1,
         tetrahedron_global_id tet2,
         tetrahedron_global_id tet3)
    : WmVol(idx, cdef, vol)
    , pNextTet()
    , pAreas()
    , pDist() {
    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    pNextTris.resize(4);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i = 0; i <= 3; ++i) {
        pNextTet[i] = nullptr;
        pNextTris[i] = nullptr;
    }
    pTets[0] = tet0;
    pTets[1] = tet1;
    pTets[2] = tet2;
    pTets[3] = tet3;

    pAreas[0] = a0;
    pAreas[1] = a1;
    pAreas[2] = a2;
    pAreas[3] = a3;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
    pDist[3] = d3;

    std::fill(pDiffBndDirection.begin(), pDiffBndDirection.end(), false);
    kprocs().resize(compdef()->countDiffs() + compdef()->countReacs());
}

////////////////////////////////////////////////////////////////////////////////

Tet::~Tet() = default;

////////////////////////////////////////////////////////////////////////////////

void Tet::checkpoint(std::fstream& cp_file) {
    WmVol::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Tet::restore(std::fstream& cp_file) {
    WmVol::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTet(uint i, Tet* t) {
    // Now adding all tets, even those from other compartments, due to the
    // diffusion boundaries
    pNextTet[i] = t;

    // if (pNextTris[i] != 0) CLOG(INFO, "general_log") << "WARNING: writing over
    // nextTri index " << i;
    pNextTris[i] = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setDiffBndDirection(uint i) {
    pDiffBndDirection.at(i) = true;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(Tri* /*t*/) {
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(uint i, Tri* t) {
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    pNextTet[i] = nullptr;
    pNextTris[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupKProcs(Tetexact* tex) {
    uint j = 0;

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (auto i: solver::reac_local_id::range(nreacs)) {
        auto& rdef = compdef()->reacdef(i);
        auto* r = new Reac(&rdef, this);
        kprocs()[j++] = r;
        tex->addKProc(r);
    }

    // Create diffusion kproc's.
    // NOTE: The order is important here- diffs should come after reacs,
    // because diffs will not be stored in WmVols and the Comp will call the
    // parent method often.
    uint ndiffs = compdef()->countDiffs();
    for (auto i: solver::diff_local_id::range(ndiffs)) {
        auto& ddef = compdef()->diffdef(i);
        auto* d = new Diff(&ddef, this);
        kprocs()[j++] = d;
        tex->addKProc(d);
    }
}

////////////////////////////////////////////////////////////////////////////////

Diff& Tet::diff(solver::diff_local_id lidx) const {
    AssertLog(lidx < compdef()->countDiffs());
    return *dynamic_cast<Diff*>(pKProcs[compdef()->countReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

int Tet::getTetDirection(tetrahedron_global_id tidx) const {
    for (uint i = 0; i < 4; i++) {
        if (pTets[i] == tidx) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

}  // namespace steps::tetexact
