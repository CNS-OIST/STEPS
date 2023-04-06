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


// Standard library & STL headers.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <functional>
#include <iostream>
#include <iostream>
#include <limits>
#include <sstream>
#include <vector>

// STEPS headers.
#include "tet.hpp"
#include "diff.hpp"
#include "reac.hpp"
#include "tetopsplit.hpp"
#include "tri.hpp"
#include "solver/diffdef.hpp"
#include "solver/reacdef.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::Tet::Tet
  (
    tetrahedron_id_t idx, solver::Compdef *cdef, double vol,
    double a0, double a1, double a2, double a3,
    double d0, double d1, double d2, double d3,
    tetrahedron_id_t tet0, tetrahedron_id_t tet1, tetrahedron_id_t tet2, tetrahedron_id_t tet3, int rank, int host_rank
  )
: WmVol(idx, cdef, vol, rank, host_rank)
, pTets({{tet0, tet1, tet2, tet3}})
, pNextTet()
, pAreas()
, pDist()
, pPoolOccupancy(nullptr)
, pLastUpdate(nullptr)
{
    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    pNextTris.resize(4);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i=0; i <= 3; ++i)
    {
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
    pPoolOccupancy = new double[nspecs];
    std::fill_n(pPoolOccupancy, nspecs, 0.0);
    pLastUpdate = new double[nspecs];
    std::fill_n(pLastUpdate, nspecs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

smtos::Tet::~Tet()
{
    delete[] pPoolOccupancy;
    delete[] pLastUpdate;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(pDiffBndDirection), sizeof(bool) * 4);
    WmVol::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(pDiffBndDirection), sizeof(bool) * 4);
    WmVol::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setNextTet(uint i, smtos::Tet * t)
{
    // Now adding all tets, even those from other compartments, due to the diffusion boundaries
    pNextTet[i] = t;

    //if (pNextTris[i] != 0) CLOG(INFO, "general_log") << "WARNING: writing over nextTri index " << i;
    pNextTris[i] = 0;

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setDiffBndDirection(uint i)
{
    AssertLog(i < 4);

    pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setNextTri(smtos::Tri */*t*/)
{
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setNextTri(uint i, smtos::Tri * t)
{
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    // This is too common now to include this message- for any internal patch this happens
    //if (pNextTet[i] != 0) CLOG(INFO, "general_log") << "WARNING: writing over nextTet index " << i;

    pNextTet[i] = nullptr;
    pNextTris[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setupKProcs(smtos::TetOpSplitP * tex)
{
    startKProcIdx = tex->countKProcs();
    uint j = 0;
    uint nreacs = compdef()->countReacs();
    uint ndiffs = compdef()->countDiffs();
    nKProcs = nreacs + ndiffs;
    // if in host create KProc
    if (hostRank == myRank) {
        // Create reaction kproc's.

        pKProcs.resize(nreacs + ndiffs);
        for (uint i = 0; i < nreacs; ++i)
        {
            ssolver::Reacdef * rdef = compdef()->reacdef(i);
            smtos::Reac * r = new smtos::Reac(rdef, this);
            pKProcs[j++] = r;
            uint idx = tex->addKProc(r);
            r->setSchedIDX(idx);
        }

        for (uint i = 0; i < ndiffs; ++i)
        {
            ssolver::Diffdef * ddef = compdef()->diffdef(i);
            auto * d = new smtos::Diff(ddef, this);
            kprocs()[j++] = d;
            uint idx = tex->addKProc(d);
            d->setSchedIDX(idx);
            tex->addDiff(d);
        }
    }
    // else just record the idx
    else {
        pKProcs.resize(0);

        for (uint i = 0; i < nKProcs; ++i)
        {
            tex->addKProc(nullptr);
        }
    }

    // Create diffusion kproc's.
    // NOTE: The order is important here- diffs should come after reacs,
    // because diffs will not be stored in WmVols and the Comp will call the
    // parent method often.

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setupDeps()
{
    super_type::setupDeps();

    bool has_remote_neighbors = false;
    const uint nspecs = compdef()->countSpecs();
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        smtos::Tet * next = nextTet(i);
        if (next == nullptr) { continue;
}
        if (next->getHost() != getHost()) {
            has_remote_neighbors = true;
            break;
        }
    }

    if (has_remote_neighbors == false) {
        localSpecUpdKProcs.clear();
        return;
    }

    localSpecUpdKProcs.resize(nspecs);
    for (uint slidx = 0; slidx < nspecs; slidx++) {
        uint sgidx = compdef()->specL2G(slidx);
        // search dependency for kprocs in this tet
        uint nkprocs = countKProcs();

        for (uint k = 0; k < nkprocs; k++)
        {
            if (KProcDepSpecTet(k, this, sgidx)) {
                localSpecUpdKProcs[slidx].push_back(getKProc(k));
            }
        }

        // search dependency for kprocs in neighboring tris
        for (uint i = 0; i < 4; ++i)
        {
            smtos::Tri * next = nextTri(i);
            if (next == nullptr) { continue;
}

            // next tri has to be in the same host to prevent
            // cross process surface reaction
            if (next->getHost() != getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next->idx() << " and its compartment tetrahedron " << idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = next->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++)
            {
                if (next->KProcDepSpecTet(sk, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(next->getKProc(sk));
                }
            }
        }
    }

}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Tet::KProcDepSpecTet(uint kp_lidx, smtos::WmVol* kp_container,  uint spec_gidx)
{
    // if kp is reaction
    uint remain = kp_lidx;
    if (remain < compdef()->countReacs()) {
        if (kp_container != this) return false;
        ssolver::Reacdef * rdef = compdef()->reacdef(remain);
        return rdef->dep(spec_gidx) != 0;
    }
    remain -= compdef()->countReacs();
    AssertLog(remain < compdef()->countDiffs());
    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        if (kp_container != this) return false;
        return spec_gidx == compdef()->diffdef(remain)->lig();
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Tet::KProcDepSpecTri(uint /*kp_lidx*/, smtos::Tri* /*kp_container*/, uint /*spec_gidx*/)
{
    // Reac and Diff never depend on species on triangle
    return false;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Diff * smtos::Tet::diff(uint lidx) const
{
    AssertLog(lidx < compdef()->countDiffs());
    return dynamic_cast<smtos::Diff*>(pKProcs[compdef()->countReacs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

 int smtos::Tet::getTetDirection(tetrahedron_id_t tidx)
{
    for (uint i = 0; i < 4; i++) {
        if (pTets[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setCount(uint lidx, uint count, double period)
{
	// Count has changed, need to correct pool factor


    // This function only apply global molecule change, i.e., called by sim.setTetCount,
	// thus skips _applyRemoteMoleculeChanges routing,
    // and only register itself if it requires sync to sreac tris.
    AssertLog(lidx < compdef()->countSpecs());
    uint oldcount = pPoolCount[lidx];
	pPoolCount[lidx] = count;

	if (period == 0.0) { return;
}

	// Count has changed,
	double lastupdate = pLastUpdate[lidx];
	AssertLog(period >= lastupdate);
	pPoolOccupancy[lidx] += oldcount*(period-lastupdate);

	pLastUpdate[lidx] = period;

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::incCount(uint lidx, int inc, double period, bool local_change)
{
    AssertLog(lidx < compdef()->countSpecs());


    // remote change caused by diffusion
    if (hostRank != myRank && !local_change)
    {
        if (inc <= 0) {
            std::ostringstream os;
            os << "Try to change molecule " << lidx << " by " << inc << "\n";
            os << "Fail because molecule change of receiving end should always be non-negative.\n";
            ProgErrLog(os.str());
        }

        bufferLocations[lidx] = pSol->registerRemoteMoleculeChange(hostRank, bufferLocations[lidx], SUB_TET, pIdx.get(), lidx, inc);
        // does not need to check sync
    }
    // local change
    else {
        double oldcount = pPoolCount[lidx];
		AssertLog(oldcount + inc >= 0.0);
		pPoolCount[lidx] += inc;


		if (period == 0.0 || local_change) { return;
}
		// Count has changed,
		double lastupdate = pLastUpdate[lidx];
		AssertLog(period >= lastupdate);
		pPoolOccupancy[lidx] += oldcount*(period-lastupdate);

		pLastUpdate[lidx] = period;



    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::resetPoolOccupancy()
{
    uint nspecs = compdef()->countSpecs();
    std::fill_n(pPoolOccupancy, nspecs, 0.0);
    std::fill_n(pLastUpdate, nspecs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::Tet::getSpecUpdKProcs(uint slidx)
{
    return localSpecUpdKProcs[slidx];
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::repartition(smtos::TetOpSplitP * tex, int rank, int host_rank)
{
    myRank = rank;
    hostRank = host_rank;

    // Delete reaction rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;

    setupKProcs(tex);
    localSpecUpdKProcs.clear();
    bufferLocations.clear();
}
////////////////////////////////////////////////////////////////////////////////

void smtos::Tet::setupBufferLocations()
{
    uint nspecs = pCompdef->countSpecs();
    bufferLocations.assign(nspecs, std::numeric_limits<uint>::max());
}

////////////////////////////////////////////////////////////////////////////////
// END
