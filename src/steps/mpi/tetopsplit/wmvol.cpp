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

#include "wmvol.hpp"

// STEPS headers.
#include "diff.hpp"
#include "reac.hpp"
#include "tetopsplit.hpp"
#include "tri.hpp"

#include "math/constants.hpp"
#include "solver/compdef.hpp"
#include "solver/diffdef.hpp"
#include "solver/reacdef.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::mpi::tetopsplit {


////////////////////////////////////////////////////////////////////////////////

WmVol::WmVol(tetrahedron_global_id idx, solver::Compdef* cdef, double vol, int rank, int host_rank)
    : pIdx(idx)
    , pCompdef(cdef)
    , pVol(vol)
    , myRank(rank)
    , hostRank(host_rank) {
    AssertLog(pCompdef != nullptr);
    AssertLog(pVol > 0.0);

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount.container().resize(nspecs);
    pPoolFlags.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

WmVol::~WmVol() {
    // Delete reaction rules.
    auto e = pKProcs.end();
    for (auto i = pKProcs.begin(); i != e; ++i) {
        delete *i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolFlags);
    util::restore(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setNextTri(Tri* t) {
    pNextTris.push_back(t);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setupKProcs(TetOpSplitP* tex) {
    startKProcIdx = tex->countKProcs();
    uint j = 0;
    nKProcs = compdef()->countReacs();
    // if in host create KProc
    if (hostRank == myRank) {
        // Create reaction kproc's.
        pKProcs.resize(nKProcs);
        for (auto i: solver::reac_local_id::range(nKProcs)) {
            auto& rdef = compdef()->reacdef(i);
            auto* r = new Reac(&rdef, this);
            pKProcs[j++] = r;
            solver::kproc_global_id idx = tex->addKProc(r);
            r->setSchedIDX(idx);
        }
    }
    // else just record the idx
    else {
        pKProcs.resize(0);

        for (uint i = 0; i < nKProcs; ++i) {
            tex->addKProc(nullptr);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setupDeps() {
    if (myRank != hostRank) {
        return;
    }
    for (auto& kp: pKProcs) {
        kp->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////

bool WmVol::KProcDepSpecTet(uint kp_lidx, WmVol* kp_container, solver::spec_global_id spec_gidx) {
    // for wmv it is always reaction so no need to check kproc type
    if (kp_container != this) {
        return false;
    }
    const auto& rdef = compdef()->reacdef(solver::reac_local_id(kp_lidx));
    return rdef.dep(spec_gidx) != 0;
}

////////////////////////////////////////////////////////////////////////////////

bool WmVol::KProcDepSpecTri(uint /*kp_lidx*/,
                            Tri* /*kp_container*/,
                            solver::spec_global_id /*spec_gidx*/) {
    // Reac never depends on species on triangle
    return false;
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

void WmVol::setCount(solver::spec_local_id lidx, uint count, double /*period*/) {
    pPoolCount.at(lidx) = count;
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::incCount(solver::spec_local_id lidx, int inc, double /*period*/, bool local_change) {
    // remote change
    if (hostRank != myRank && !local_change) {
        std::ostringstream os;
        os << "Remote WmVol update is not implemented.\n";
        NotImplErrLog(os.str());
    }
    // local change
    else {
        double oldcount = pPoolCount.at(lidx);
        AssertLog(oldcount + inc >= 0.0);
        pPoolCount[lidx] += inc;
    }
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
    AssertLog(lidx < compdef()->countReacs());
    return *dynamic_cast<Reac*>(pKProcs[lidx.get()]);
}
////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool WmVol::getInHost() const {
    return hostRank == myRank;
}

////////////////////////////////////////////////////////////////////////////////
void WmVol::setHost(int host, int rank) {
    hostRank = host;
    myRank = rank;
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::setSolver(mpi::tetopsplit::TetOpSplitP* solver) {
    pSol = solver;
}

////////////////////////////////////////////////////////////////////////////////

TetOpSplitP* WmVol::solver() const {
    return pSol;
}

////////////////////////////////////////////////////////////////////////////////

double WmVol::getPoolOccupancy(solver::spec_local_id /*lidx*/) const {
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double WmVol::getLastUpdate(solver::spec_local_id /*lidx*/) const {
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::resetPoolOccupancy() {
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void WmVol::repartition(TetOpSplitP* tex, int rank, int host_rank) {
    myRank = rank;
    hostRank = host_rank;

    // Delete reaction rules.
    auto e = pKProcs.end();
    for (auto i = pKProcs.begin(); i != e; ++i) {
        delete *i;
    }

    setupKProcs(tex);
}

}  // namespace steps::mpi::tetopsplit
