/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

#include <mpi.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/reacdef.hpp"
#include "steps/mpi/mpi_common.hpp"
#include "steps/mpi/tetopsplit/diff.hpp"
#include "steps/mpi/tetopsplit/reac.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
#include "steps/mpi/tetopsplit/wmvol.hpp"

#include "third_party/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::WmVol::WmVol
(
    uint idx, solver::Compdef * cdef, double vol, int rank, int host_rank
)
: pIdx(idx)
, pCompdef(cdef)
, pVol(vol)
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
, pNextTris()
, myRank(rank)
, hostRank(host_rank)
{
    assert(pCompdef != 0);
    assert (pVol > 0.0);

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);

}

////////////////////////////////////////////////////////////////////////////////

smtos::WmVol::~WmVol(void)
{
    // Delete species pool information.
    delete[] pPoolCount;
    delete[] pPoolFlags;

    // Delete reaction rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::checkpoint(std::fstream & cp_file)
{
    uint nspecs = compdef()->countSpecs();
    cp_file.write((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.write((char*)pPoolFlags, sizeof(uint) * nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::restore(std::fstream & cp_file)
{
    uint nspecs = compdef()->countSpecs();
    cp_file.read((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.read((char*)pPoolFlags, sizeof(uint) * nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setNextTri(smtos::Tri * t)
{
    uint index = pNextTris.size();
    pNextTris.push_back(t);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setupKProcs(smtos::TetOpSplitP * tex)
{
    startKProcIdx = tex->countKProcs();
    uint j = 0;
    nKProcs = compdef()->countReacs();
    // if in host create KProc
    if (hostRank == myRank) {
        // Create reaction kproc's.
        pKProcs.resize(nKProcs);
        for (uint i = 0; i < nKProcs; ++i)
        {
            ssolver::Reacdef * rdef = compdef()->reacdef(i);
            smtos::Reac * r = new smtos::Reac(rdef, this);
            pKProcs[j++] = r;
            uint idx = tex->addKProc(r, hostRank);
            r->setSchedIDX(idx);
        }
    }
    // else just record the idx
    else {
        pKProcs.resize(0);
        
        for (uint i = 0; i < nKProcs; ++i)
        {
            uint idx = tex->addKProc(NULL, hostRank);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setupDeps(void)
{
    if (myRank != hostRank) return;
    for (auto kp : pKProcs) {
        kp->setupDeps();
    }
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::WmVol::KProcDepSpecTet(uint kp_lidx, smtos::WmVol* kp_container,  uint spec_gidx)
{
    // for wmv it is always reaction so no need to check kproc type
    if (kp_container != this) return false;
    ssolver::Reacdef * rdef = compdef()->reacdef(kp_lidx);
    return rdef->dep(spec_gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::WmVol::KProcDepSpecTri(uint kp_lidx, smtos::Tri* kp_container, uint spec_gidx)
{
    // Reac never depends on species on triangle
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::reset(void)
{
    uint nspecs = compdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(),
            std::mem_fun(&smtos::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

double smtos::WmVol::conc(uint gidx) const
{
    uint lspidx = compdef()->specG2L(gidx);
    double n = pPoolCount[lspidx];
    return (n/(1.0e3*pVol*steps::math::AVOGADRO));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setCount(uint lidx, uint count, double period)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Set spec (local id): " << lidx << " to " << count << "\n";
    #endif

    assert (lidx < compdef()->countSpecs());
    pPoolCount[lidx] = count;

}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::incCount(uint lidx, int inc, double period)
{
    #ifdef MPI_DEBUG
    //CLOG(DEBUG, "mpi_debug") << "Change spec (local id): " << lidx << " by " << inc << "\n";
    #endif

    assert (lidx < compdef()->countSpecs());
    
    
    // remote change
    if (hostRank != myRank)
    {
        std::ostringstream os;
        os << "Remote WmVol update is not implemented.\n";
        throw steps::NotImplErr(os.str());
    }
    // local change
    else {
        double oldcount = pPoolCount[lidx];
		assert(oldcount + inc >= 0.0);
		pPoolCount[lidx] += inc;
		
		#ifdef MPI_DEBUG
        //CLOG(DEBUG, "mpi_debug") << "local count: " << pPoolCount[lidx] << "\n";
        #endif
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

smtos::Reac * smtos::WmVol::reac(uint lidx) const
{
    assert(lidx < compdef()->countReacs());
    return dynamic_cast<smtos::Reac*>(pKProcs[lidx]);
}
////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool smtos::WmVol::getInHost(void)
{
    return (hostRank == myRank);
}

////////////////////////////////////////////////////////////////////////////////
void smtos::WmVol::setHost(int host, int rank)
{
    hostRank = host;
    myRank = rank;
}

////////////////////////////////////////////////////////////////////////////////
/*
void smtos::WmVol::addSyncHost(int host)
{
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "add sync host " << host << ".\n";
    #endif
    syncHosts.insert(host);
}

*/
////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::setSolver(steps::mpi::tetopsplit::TetOpSplitP* solver)
{
    pSol = solver;
}

////////////////////////////////////////////////////////////////////////////////

smtos::TetOpSplitP* smtos::WmVol::solver(void)
{
    return pSol;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::WmVol::getPoolOccupancy(uint lidx)
{
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::WmVol::getLastUpdate(uint lidx)
{
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::resetPoolOccupancy(void)
{
    assert (false);
}

////////////////////////////////////////////////////////////////////////////////


void smtos::WmVol::repartition(smtos::TetOpSplitP * tex, int rank, int host_rank)
{
    myRank = rank;
    hostRank = host_rank;
    
    // Delete reaction rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
    
    setupKProcs(tex);
}

////////////////////////////////////////////////////////////////////////////////
/*
void smtos::WmVol::sendSyncPools(void)
{
    uint nspecs = pCompdef->countSpecs();
    smtos::Tet* upcast = dynamic_cast<smtos::Tet*> (this);
    uint info[2];
    info[0] = smtos::SUB_WM;
    if (upcast != NULL) info[0] = smtos::SUB_TET;
    info[1] = pIdx;
    for (auto & host : syncHosts) {
        #ifdef MPI_DEBUG
        CLOG(DEBUG, "mpi_debug") << "Send pool " << pIdx << " to " << host << ".\n";
		std::vector<uint> temp;
	    for(uint s = 0; s < nspecs; s++) {
	        temp.push_back(pPoolCount[s]);
	    }
        CLOG(DEBUG, "mpi_debug") << temp << ".\n";
        #endif
        MPI_Send(info, 2, MPI_UNSIGNED, host, OPSPLIT_COUNT_SYNC_INFO, MPI_COMM_WORLD);
        MPI_Send(pPoolCount, nspecs, MPI_UNSIGNED, host, OPSPLIT_COUNT_SYNC_DATA, MPI_COMM_WORLD);
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::WmVol::recvSyncPools(int source)
{
    #ifdef MPI_DEBUG
    CLOG(DEBUG, "mpi_debug") << "Recv pool " << pIdx << " from " << source << ".\n";
    #endif
    uint nspecs = pCompdef->countSpecs();
    MPI_Recv(pPoolCount, nspecs, MPI_UNSIGNED, source, OPSPLIT_COUNT_SYNC_DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    #ifdef MPI_DEBUG
    std::vector<uint> temp;
	for(uint s = 0; s < nspecs; s++) {
	    temp.push_back(pPoolCount[s]);
	}
    CLOG(DEBUG, "mpi_debug") << temp << ".\n";
    #endif
}
*/

/////
// MPISTEPS
////////////////////////////////////////////////////////////////////////////////

// END
