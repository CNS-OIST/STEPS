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
#include "tet.hpp"

#include "diff.hpp"
#include "reac.hpp"
#include "solver/diffdef.hpp"
#include "solver/reacdef.hpp"
#include "tetopsplit.hpp"
#include "tri.hpp"

// logging
#include "util/error.hpp"

#include "util/checkpointing.hpp"

namespace steps::mpi::tetopsplit {


////////////////////////////////////////////////////////////////////////////////

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
         tetrahedron_global_id tet3,
         int rank,
         int host_rank)
    : WmVol(idx, cdef, vol, rank, host_rank)
    , pTets({{tet0, tet1, tet2, tet3}})
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

    std::fill_n(pDiffBndDirection, 4, false);

    uint nspecs = compdef()->countSpecs();
    pPoolOccupancy.container().resize(nspecs);
    pLastUpdate.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void Tet::checkpoint(std::fstream& cp_file) {
    WmVol::checkpoint(cp_file);
    // NOTE not checkpointing pPoolOccupancy or pLastUpdate because these should be reset at time of
    // call to checkpoint
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
    AssertLog(i < 4);

    pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(Tri* /*t*/) {
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTri(uint i, Tri* t) {
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    // This is too common now to include this message- for any internal patch this
    // happens
    // if (pNextTet[i] != 0) CLOG(INFO, "general_log") << "WARNING: writing over
    // nextTet index " << i;

    pNextTet[i] = nullptr;
    pNextTris[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupKProcs(TetOpSplitP* tex) {
    startKProcIdx = tex->countKProcs();
    uint j = 0;
    uint nreacs = compdef()->countReacs();
    uint ndiffs = compdef()->countDiffs();
    nKProcs = nreacs + ndiffs;
    // if in host create KProc
    if (hostRank == myRank) {
        // Create reaction kproc's.

        pKProcs.resize(nreacs + ndiffs);
        for (auto i: solver::reac_local_id::range(nreacs)) {
            auto& rdef = compdef()->reacdef(i);
            Reac* r = new Reac(&rdef, this);
            pKProcs[j++] = r;
            solver::kproc_global_id idx = tex->addKProc(r);
            r->setSchedIDX(idx);
        }

        for (auto i: solver::diff_local_id::range(ndiffs)) {
            auto& ddef = compdef()->diffdef(i);
            auto* d = new Diff(&ddef, this);
            kprocs()[j++] = d;
            solver::kproc_global_id idx = tex->addKProc(d);
            d->setSchedIDX(idx);
            tex->addDiff(d);
        }
    }
    // else just record the idx
    else {
        pKProcs.resize(0);

        for (uint i = 0; i < nKProcs; ++i) {
            tex->addKProc(nullptr);
        }
    }

    // Create diffusion kproc's.
    // NOTE: The order is important here- diffs should come after reacs,
    // because diffs will not be stored in WmVols and the Comp will call the
    // parent method often.
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupDeps() {
    super_type::setupDeps();

    bool has_remote_neighbors = false;
    const uint nspecs = compdef()->countSpecs();
    for (uint i = 0; i < 4; ++i) {
        // Fetch next tetrahedron, if it exists.
        Tet* next = nextTet(i);
        if (next == nullptr) {
            continue;
        }
        if (next->getHost() != getHost()) {
            has_remote_neighbors = true;
            break;
        }
    }

    if (has_remote_neighbors == false) {
        localSpecUpdKProcs.container().clear();
        return;
    }

    localSpecUpdKProcs.container().resize(nspecs);
    for (auto slidx: solver::spec_local_id::range(nspecs)) {
        solver::spec_global_id sgidx = compdef()->specL2G(slidx);
        // search dependency for kprocs in this tet
        uint nkprocs = countKProcs();

        for (uint k = 0; k < nkprocs; k++) {
            if (KProcDepSpecTet(k, this, sgidx)) {
                localSpecUpdKProcs[slidx].push_back(getKProc(k));
            }
        }

        // search dependency for kprocs in neighboring tris
        for (uint i = 0; i < 4; ++i) {
            Tri* next = nextTri(i);
            if (next == nullptr) {
                continue;
            }

            // next tri has to be in the same host to prevent
            // cross process surface reaction
            if (next->getHost() != getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next->idx() << " and its compartment tetrahedron "
                   << idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = next->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                if (next->KProcDepSpecTet(sk, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(next->getKProc(sk));
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Tet::KProcDepSpecTet(uint kp_lidx, WmVol* kp_container, solver::spec_global_id spec_gidx) {
    // if kp is reaction
    uint remain = kp_lidx;
    if (remain < compdef()->countReacs()) {
        if (kp_container != this) {
            return false;
        }
        auto& rdef = compdef()->reacdef(solver::reac_local_id(remain));
        return rdef.dep(spec_gidx) != 0;
    }
    remain -= compdef()->countReacs();
    AssertLog(remain < compdef()->countDiffs());
    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        if (kp_container != this) {
            return false;
        }
        return spec_gidx == compdef()->diffdef(solver::diff_local_id(remain)).lig();
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool Tet::KProcDepSpecTri(uint /*kp_lidx*/,
                          Tri* /*kp_container*/,
                          solver::spec_global_id /*spec_gidx*/) {
    // Reac and Diff never depend on species on triangle
    return false;
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
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setCount(solver::spec_local_id lidx, uint count, double period) {
    // Count has changed, need to correct pool factor

    // This function is used for updates that do not require remote sync,
    // such as user setCount or reaction

    AssertLog(lidx < compdef()->countSpecs());
    uint oldcount = pPoolCount[lidx];
    pPoolCount[lidx] = count;

    if (period == 0.0) {
        return;
    }

    // Count has changed,
    double lastupdate = pLastUpdate[lidx];
    AssertLog(period >= lastupdate);
    pPoolOccupancy[lidx] += oldcount * (period - lastupdate);

    pLastUpdate[lidx] = period;
}

////////////////////////////////////////////////////////////////////////////////

void Tet::incCount(solver::spec_local_id lidx, int inc, double period, bool local_change) {
    AssertLog(lidx < compdef()->countSpecs());

    // remote change caused by diffusion
    if (hostRank != myRank && !local_change) {
        if (inc <= 0) {
            std::ostringstream os;
            os << "Try to change molecule " << lidx << " by " << inc << "\n";
            os << "Fail because molecule change of receiving end should always be "
                  "non-negative.\n";
            ProgErrLog(os.str());
        }

        bufferLocations[lidx] = pSol->registerRemoteMoleculeChange(
            hostRank, bufferLocations[lidx], SUB_TET, pIdx.get(), lidx, inc);
        // does not need to check sync
    }
    // local change
    else {
        double oldcount = pPoolCount[lidx];
        AssertLog(oldcount + inc >= 0.0);
        pPoolCount[lidx] += inc;

        if (period == 0.0 || local_change) {
            return;
        }
        // Count has changed,
        double lastupdate = pLastUpdate[lidx];
        AssertLog(period >= lastupdate);
        pPoolOccupancy[lidx] += oldcount * (period - lastupdate);
        pLastUpdate[lidx] = period;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tet::resetPoolOccupancy() {
    std::fill(pPoolOccupancy.begin(), pPoolOccupancy.end(), 0.0);
    std::fill(pLastUpdate.begin(), pLastUpdate.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& Tet::getSpecUpdKProcs(solver::spec_local_id slidx) {
    return localSpecUpdKProcs[slidx];
}

////////////////////////////////////////////////////////////////////////////////

void Tet::repartition(TetOpSplitP* tex, int rank, int host_rank) {
    myRank = rank;
    hostRank = host_rank;

    // Delete reaction rules.
    auto e = pKProcs.end();
    for (auto i = pKProcs.begin(); i != e; ++i) {
        delete *i;
    }

    setupKProcs(tex);
    localSpecUpdKProcs.container().clear();
    bufferLocations.container().clear();
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupBufferLocations() {
    uint nspecs = pCompdef->countSpecs();
    bufferLocations.container().assign(nspecs, std::numeric_limits<uint>::max());
}

}  // namespace steps::mpi::tetopsplit
