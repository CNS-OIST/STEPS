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
#include <vector>
#include <set>
#include <fstream>
#include <sstream>
#include <random>

#include <mpi.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/mpi/tetopsplit/diff.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
#include "third_party/easylogging++.h"

#include <iostream>
#include <mpi.h>

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

smtos::Diff::Diff(ssolver::Diffdef * ddef, smtos::Tet * tet)
: KProc()
, pDiffdef(ddef)
, pTet(tet)
, localAllUpdVec()
, remoteAllUpdVec()
, pScaledDcst(0.0)
, pDcst(0.0)
, pNonCDFSelector()
, pNeighbCompLidx()
, pEmptyvec()
, idxEmptyvec()
, pDirections()
, pNdirections(0)
{
    assert(pDiffdef != 0);
    assert(pTet != 0);
    type = KP_DIFF;
    smtos::Tet * next[4] =
    {
        pTet->nextTet(0),
        pTet->nextTet(1),
        pTet->nextTet(2),
        pTet->nextTet(3)
    };

    ligGIdx = pDiffdef->lig();
    ssolver::Compdef * cdef = pTet->compdef();
    lidxTet = cdef->specG2L(ligGIdx);

    for (uint i = 0; i < 4; ++i) { pDiffBndDirection[i] = pTet->getDiffBndDirection(i);}

    uint gidx = pDiffdef->lig();
    for (uint i = 0; i < 4; ++i)
    {
        if (next[i] == 0)
        {
            pNeighbCompLidx[i] = -1;
            continue;
        }
        else
        {
            pNeighbCompLidx[i] = next[i]->compdef()->specG2L(ligGIdx);
        }
    }

    // Precalculate part of the scaled diffusion constant.
    uint ldidx = pTet->compdef()->diffG2L(pDiffdef->gidx());
    double dcst = pTet->compdef()->dcst(ldidx);
    pDcst = dcst;

    double d[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pDiffBndDirection[i] == true)
            {
                if (pDiffBndActive[i]) d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == cdef)    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0;
            }
        }
    }
    // Compute scaled "diffusion constant".
    for (uint i = 0; i < 4; ++i)
    {
        pScaledDcst += d[i];
    }

    // Should not be negative!
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
        pNonCDFSelector[3] = 0.0;

    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;
        pNonCDFSelector[3] = d[3]/pScaledDcst;
        
        pNdirections = 0;
        pDirections.clear();
        
        if (d[0] > 0.0)
        {
			pDirections.push_back(0);
			pNdirections+=1;
       	}
		if (d[1] > 0.0)
		{
			pDirections.push_back(1);
			pNdirections+=1;
		}
		if (d[2] > 0.0)
		{
			pDirections.push_back(2);
			pNdirections+=1;
		}
		if (d[3] > 0.0)
		{
			pDirections.push_back(3);
			pNdirections+=1;
		}

		//pMulti_n = new uint [pNdirections];
		//pDirection_prob = new double [pNdirections];
	    //for (uint i = 0; i < pNdirections; ++i) pDirection_prob[i] = d[pDirections[i]]/pScaledDcst;

    }
}

////////////////////////////////////////////////////////////////////////////////

smtos::Diff::~Diff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));

    uint n_direct_dcsts = directionalDcsts.size();
    cp_file.write((char*)&n_direct_dcsts, sizeof(uint));
    for (auto& item : directionalDcsts) {
        cp_file.write((char*)&(item.first), sizeof(uint));
        cp_file.write((char*)&(item.second), sizeof(double));
    }
    
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(n_direct_dcsts != 0, DEBUG, "steps_debug") << "Stored Directional Dcst mapping: " << directionalDcsts << "\n";
    #endif
    
    cp_file.write((char*)&pScaledDcst, sizeof(double));
    cp_file.write((char*)&pDcst, sizeof(double));
    cp_file.write((char*)pDiffBndActive, sizeof(bool) * 4);
    cp_file.write((char*)pDiffBndDirection, sizeof(bool) * 4);
    cp_file.write((char*)pNeighbCompLidx, sizeof(int) * 4);

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&rExtent, sizeof(uint));
    cp_file.read((char*)&pFlags, sizeof(uint));

    uint n_direct_dcsts = 0;
    cp_file.read((char*)&n_direct_dcsts, sizeof(uint));
    for (uint i = 0; i < n_direct_dcsts; i++) {
        uint id = 0;
        double value = 0.0;
        cp_file.read((char*)&id, sizeof(uint));
        cp_file.read((char*)&value, sizeof(double));
        directionalDcsts[id] = value;
    }
    
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG_IF(n_direct_dcsts != 0, DEBUG, "steps_debug") << "Restored Directional Dcst mapping: " << directionalDcsts << "\n";
    #endif

    cp_file.read((char*)&pScaledDcst, sizeof(double));
    cp_file.read((char*)&pDcst, sizeof(double));
    cp_file.read((char*)pDiffBndActive, sizeof(bool) * 4);
    cp_file.read((char*)pDiffBndDirection, sizeof(bool) * 4);
    cp_file.read((char*)pNeighbCompLidx, sizeof(int) * 4);

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::setupDeps(void)
{
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' tetrahedron
    //   * any neighbouring triangles
    //
    // But also in the possible 'destination' tetrahedrons (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring tetrahedrons
    //   * any neighbouring triangles of these neighbouring tets
    //
    // Since there can be no diffusion between tetrahedrons blocked by
    // a triangle, there is no need to filter out duplicate dependent
    // kprocs.

    // Search for dependencies in the 'source' tetrahedron.
    assert(pTet->getInHost());
    std::set<uint> remote;
    std::set<uint> remote_all;
    
    std::set<smtos::KProc*> local;
    std::set<smtos::KProc*> local_all;

    uint nkprocs = pTet->countKProcs();
    uint startKProcIdx = pTet->getStartKProcIdx();
    
    for (uint k = 0; k < nkprocs; k++)
    {
        if (pTet->KProcDepSpecTet(k, pTet, ligGIdx) == true) {
            local.insert(pTet->getKProc(k));
        }
    }
    // Check the neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        smtos::Tri * next = pTet->nextTri(i);
        if (next == 0) continue;
        
        // next tri has to be in the same host to prevent
        // cross process surface reaction
        if (next->getHost() != pTet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << next->idx() << " and its compartment tetrahedron " << pTet->idx()  << " belong to different hosts.\n";
            throw steps::NotImplErr(os.str());
        }
    
        nkprocs = next->countKProcs();
        startKProcIdx = next->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            if (next->KProcDepSpecTet(sk, pTet, ligGIdx) == true) {
                local.insert(next->getKProc(sk));
            }
        }
    }

    // Search for dependencies in neighbouring tetrahedra.
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        smtos::Tet * next = pTet->nextTet(i);
        if (next == 0) continue;
        if (pTet->nextTri(i) != 0) continue;

        // Copy local dependencies.
        std::set<KProc*> local2(local.begin(), local.end());
        std::set<uint> remote2;

        // Find the ones 'locally' in the next tet.
        nkprocs = next->countKProcs();
        startKProcIdx = next->getStartKProcIdx();
        
        if (next->getHost() != pTet->getHost()) {
            pTet->solver()->addNeighHost(next->getHost());
            pTet->solver()->registerBoundaryTet(next);
        }
        
        for (uint k = 0; k < nkprocs; k++)
        {
            if (next->KProcDepSpecTet(k, next, ligGIdx) == true) {
                // if next tet has different host with pTet, store dependent kp as pointer
                if (next->getHost() == pTet->getHost()) {
                    local2.insert(next->getKProc(k));
                }
                // if not store as index
                else {
                    remote2.insert(startKProcIdx + k);
                }
            }
        }

        // Find deps in neighbouring triangles in the next tet.
        // As said before, this cannot logically include the shared
        // triangle.
        for (uint j = 0; j < 4; ++j)
        {
            // Fetch next triangle, if it exists.
            smtos::Tri * next2 = next->nextTri(j);
            if (next2 == 0) continue;

            if (next2->getHost() != next->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next2->idx() << " and its compartment tetrahedron " << next->idx()  << " belong to different hosts.\n";
                throw steps::NotImplErr(os.str());
            }
            
            // Find deps.
            nkprocs = next2->countKProcs();
            startKProcIdx = next2->getStartKProcIdx();
            
            for (uint sk = 0; sk < nkprocs; sk++)
            {
                if (next2->KProcDepSpecTet(sk, next, ligGIdx) == true) {
                    if (next2->getHost() == pTet->getHost()) {
                        local2.insert(next2->getKProc(sk));
                    }
                    else {
                        remote2.insert(startKProcIdx + sk);
                    }
                }
            }
        }

        // Copy the set to the update vector.
        localUpdVec[i].assign(local2.begin(), local2.end());
        remoteUpdVec[i].assign(remote2.begin(), remote2.end());
        local_all.insert(local2.begin(), local2.end());
        remote_all.insert(remote2.begin(), remote2.end());
    }
    localAllUpdVec.assign(local_all.begin(), local_all.end());
    remoteAllUpdVec.assign(remote_all.begin(), remote_all.end());
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Diff::depSpecTet(uint gidx, smtos::WmVol * tet)
{
    if (pTet != tet) return false;
    if (gidx != ligGIdx) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Diff::depSpecTri(uint gidx, smtos::Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::reset(void)
{
    resetExtent();

    // NOTE: These must become the dcst calculation for obvious reasons
    pDiffBndActive[0] = false;
    pDiffBndActive[1] = false;
    pDiffBndActive[2] = false;
    pDiffBndActive[3] = false;

    uint ldidx = pTet->compdef()->diffG2L(pDiffdef->gidx());
    double dcst = pTet->compdef()->dcst(ldidx);
    
    // directional dcst will also be clear by setDcst
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    // pos cannot be reset because their positions in pDiffs will be used to swap
    //crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::setDiffBndActive(uint i, bool active)
{
    assert (i < 4);
    assert(pDiffBndDirection[i] == true);

    // Only need to update if the flags are changing
    if (pDiffBndActive[i] != active)
    {
        pDiffBndActive[i] = active;
        setDcst(pDcst);
    }

}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Diff::getDiffBndActive(uint i) const
{
    assert (i < 4);
    assert(pDiffBndDirection[i] == true);

    return pDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Diff::dcst(int direction)
{
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    }
    else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::setDcst(double dcst)
{
    assert(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();
    
    smtos::Tet * next[4] =
    {
    pTet->nextTet(0),
    pTet->nextTet(1),
    pTet->nextTet(2),
    pTet->nextTet(3)
    };
    
	// Reset this stuff- may have been created before, may not have been (if original dcst was 0)
    pNdirections=0;
    pDirections.clear();
    
    double d[4] = { 0.0, 0.0, 0.0, 0.0 };

    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pDiffBndDirection[i] == true)
            {
                if (pDiffBndActive[i]) d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == pTet->compdef())    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0;
            }
         }
    }

    // Compute scaled "diffusion constant".
    pScaledDcst = 0.0;

    for (uint i = 0; i < 4; ++i)
    {
        pScaledDcst += d[i];
    }
    // Should not be negative!
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
        pNonCDFSelector[3] = 0.0;
    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;
        pNonCDFSelector[3] = d[3]/pScaledDcst;

        if (d[0] > 0.0)
        {
			pDirections.push_back(0);
			pNdirections+=1;
       	}
		if (d[1] > 0.0)
		{
			pDirections.push_back(1);
			pNdirections+=1;
		}
		if (d[2] > 0.0)
		{
			pDirections.push_back(2);
			pNdirections+=1;
		}
		if (d[3] > 0.0)
		{
			pDirections.push_back(3);
			pNdirections+=1;
		}

		//pMulti_n = new uint [pNdirections];
		//pDirection_prob = new double [pNdirections];
	    //for (uint i = 0; i < pNdirections; ++i) pDirection_prob[i] = d[pDirections[i]]/pScaledDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Diff::setDirectionDcst(int direction, double dcst)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG(DEBUG, "steps_debug") << "Direction: " << direction << " dcst: " << dcst << "\n";
    #endif
    
    assert(direction < 4);
    assert(direction >= 0);
    assert(dcst >= 0.0);
    directionalDcsts[direction] = dcst;
    
    // Automatically activate boundary diffusion if necessary
    if (pDiffBndDirection[direction] == true) pDiffBndActive[direction] = true;

    smtos::Tet * next[4] =
    {
        pTet->nextTet(0),
        pTet->nextTet(1),
        pTet->nextTet(2),
        pTet->nextTet(3)
    };

    double d[4] = { 0.0, 0.0, 0.0, 0.0 };
    
    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pDiffBndDirection[i] == true)
            {
                #ifdef DIRECTIONAL_DCST_DEBUG
                CLOG(DEBUG, "steps_debug") << "enable boundary\n";
                CLOG(DEBUG, "steps_debug") << "activation: " << pDiffBndActive[i] << "\n";
                #endif
                
                if (pDiffBndActive[i]) {
                    if (directionalDcsts.find(i) != directionalDcsts.end()) {
                    
                        #ifdef DIRECTIONAL_DCST_DEBUG
                        CLOG(DEBUG, "steps_debug") << "enable boundary, use stored directional dcst: " << directionalDcsts[i] << " to compute d[" << i << "]\n";
                        #endif
                        
                        d[i] = (pTet->area(i) * directionalDcsts[i]) / (pTet->vol() * dist);
                    }
                    else
                    {
                        #ifdef DIRECTIONAL_DCST_DEBUG
                        CLOG(DEBUG, "steps_debug") << "enable boundary, use default dcst: " << pDcst << " to compute d[" << i << "]\n";
                        #endif
                        
                        d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
                    }
                }
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == pTet->compdef())    {
                    if (directionalDcsts.find(i) != directionalDcsts.end()) {
                    
                        #ifdef DIRECTIONAL_DCST_DEBUG
                        CLOG(DEBUG, "steps_debug")  << "no boundary, use stored directional dcst: " << directionalDcsts[i] << " to compute d[" << i << "]\n";
                        #endif
                        
                        d[i] = (pTet->area(i) * directionalDcsts[i]) / (pTet->vol() * dist);
                    }
                    else {
                    
                        #ifdef DIRECTIONAL_DCST_DEBUG
                        CLOG(DEBUG, "steps_debug") << "no boundary, use default dcst: " << pDcst << " to compute d[" << i << "]\n";
                        #endif
                        
                        d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
                    }
                }
                else d[i] = 0;
            }
        }
    }
    
    // Compute scaled "diffusion constant".
    pScaledDcst = 0.0;
    
    for (uint i = 0; i < 4; ++i)
    {
        pScaledDcst += d[i];
    }
    // Should not be negative!
    assert(pScaledDcst >= 0);
    // Reset this stuff- may have been created before, may not have been (if original dcst was 0)
    pNdirections=0;
    pDirections.clear();
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
        pNonCDFSelector[3] = 0.0;
    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;
        pNonCDFSelector[3] = d[3]/pScaledDcst;
        
        if (d[0] > 0.0)
        {
            pDirections.push_back(0);
            pNdirections+=1;
        }
        if (d[1] > 0.0)
        {
            pDirections.push_back(1);
            pNdirections+=1;
        }
        if (d[2] > 0.0)
        {
            pDirections.push_back(2);
            pNdirections+=1;
        }
        if (d[3] > 0.0)
        {
            pDirections.push_back(3);
            pNdirections+=1;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Diff::rate(steps::mpi::tetopsplit::TetOpSplitP * solver)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTet->pools()[lidxTet]);
    assert(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

//double smtos::Diff::getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver)
//{
//    return pScaledDcst;
//}



////////////////////////////////////////////////////////////////////////////////

int smtos::Diff::apply(steps::rng::RNG * rng)
{

	// Apply local change.
	uint * local = pTet->pools() + lidxTet;
	bool clamped = pTet->clamped(lidxTet);

	if (clamped == false)
	{
		//assert(*local > 0);
		if (*local == 0)
		{
			return -2;
		}
	}

	// Apply change in next voxel: select a direction.
	double sel = rng->getUnfEE();

	if (sel < pNonCDFSelector[0])
	{

		// Direction 1.
		smtos::Tet * nexttet = pTet->nextTet(0);
		// If there is no next tet 0, pCDFSelector[0] should be zero
		// So we can assert that nextet 0 does indeed exist
		assert (nexttet != 0);
		assert(pNeighbCompLidx[0] > -1);

		if (nexttet->clamped(pNeighbCompLidx[0]) == false)
		{
			nexttet->incCount(pNeighbCompLidx[0],1);
		}
		if (clamped == false) {pTet->incCount(lidxTet, -1); }

		rExtent++;

		return 0;
	}
	else if (sel < pNonCDFSelector[0] + pNonCDFSelector[1])
	{

		// Direction 2.
		smtos::Tet * nexttet = pTet->nextTet(1);
		// If there is no next tet 1, pCDFSelector[1] should be zero
		// So we can assert that nextet 1 does indeed exist
		assert (nexttet != 0);
		assert(pNeighbCompLidx[1] > -1);

		if (nexttet->clamped(pNeighbCompLidx[1]) == false)
		{
			nexttet->incCount(pNeighbCompLidx[1],1);
		}

		if (clamped == false) {pTet->incCount(lidxTet, -1); }

		rExtent++;

		return 1;
	}
	else if (sel < pNonCDFSelector[0] + pNonCDFSelector[1] + pNonCDFSelector[2])
	{

		// Direction 3.
		smtos::Tet * nexttet = pTet->nextTet(2);
		// If there is no next tet 2, pCDFSelector[2] should be zero
		// So we can assert that nextet 2 does indeed exist
		assert (nexttet != 0);
		assert(pNeighbCompLidx[2] > -1);
		if (nexttet->clamped(pNeighbCompLidx[2]) == false)
		{
			nexttet->incCount(pNeighbCompLidx[2],1);
		}

		if (clamped == false) {pTet->incCount(lidxTet, -1); }

		rExtent++;
		return 2;
	}
	else
	{

		// Direction 4.
		smtos::Tet * nexttet = pTet->nextTet(3);
		// If there is no next tet 3, pCDFSelector[3] should be zero
		// So we can assert that nextet 3 does indeed exist
		assert (nexttet != 0);
		assert(pNeighbCompLidx[3] > -1);

		if (nexttet->clamped(pNeighbCompLidx[3]) == false)
		{
			nexttet->incCount(pNeighbCompLidx[3],1);
		}
		if (clamped == false) {pTet->incCount(lidxTet, -1); }

		rExtent++;

		return 3;

	}

	// This should never happen!
	std::cerr << "Cannot find a suitable direction for diffusion!\n";
	throw;
	return -1;
}

///////////////////////////////////////////////////////////////////////////////

int smtos::Diff::apply(steps::rng::RNG * rng, uint nmolcs)
{
	// Apply local change.
	uint * local = pTet->pools() + lidxTet;
	bool clamped = pTet->clamped(lidxTet);

	if (clamped == false)
	{
		//assert(*local > 0);
        if (*local == 0)
        {
            return -2;
        }
    }

	assert(pNdirections >= 1);

    // Multinomial by stl
	uint molcs_moved = 0;
    for (uint i=0; i< pNdirections-1; ++i)
    {
		uint direction = pDirections[i];
        double chance=0.0;

        double sump = 0.0;
        for (uint j=0; j<i; ++j) sump +=pNonCDFSelector[pDirections[j]];

        chance = pNonCDFSelector[direction]/(1.0-sump);

        unsigned int max_molcs= nmolcs-molcs_moved;

        // This is really important: if chance is a fraction over 1.0 the binomial
        // cycles to the largest unsigned int
        if (chance >= 1.0) chance=1.0;
        
        uint molcsthisdir = rng->getBinom(max_molcs, chance);
        
        if (molcsthisdir)
        {
        	smtos::Tet * nexttet = pTet->nextTet(direction);

        	assert (nexttet != 0);
        	assert(pNeighbCompLidx[direction] > -1);

        	if (nexttet->clamped(pNeighbCompLidx[direction]) == false)
        	{
        		nexttet->incCount(pNeighbCompLidx[direction],molcsthisdir);
        	}

        	molcs_moved+=molcsthisdir;
        }
        if (molcs_moved == nmolcs) break;

    }
    //last direction, chance =1
    uint direction = pDirections[pNdirections-1];
    int molcsthisdir = nmolcs-molcs_moved;
    if (molcsthisdir)
    {
        	smtos::Tet * nexttet = pTet->nextTet(direction);

        	assert (nexttet != 0);
        	assert(pNeighbCompLidx[direction] > -1);

        	if (nexttet->clamped(pNeighbCompLidx[direction]) == false)
        	{
        		nexttet->incCount(pNeighbCompLidx[direction],molcsthisdir);
        	}

        	molcs_moved+=molcsthisdir;
    }

	assert(molcs_moved == nmolcs);

	if (clamped == false) {pTet->incCount(lidxTet, -nmolcs); }

	rExtent+=nmolcs;
	return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::Diff::getRemoteUpdVec(int direction)
{
    if (direction == -1) return remoteAllUpdVec;
    else if (direction == -2) return idxEmptyvec;
    else return remoteUpdVec[direction];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::Diff::getLocalUpdVec(int direction)
{
    if (direction == -1) return localAllUpdVec;
    else if (direction == -2) return pEmptyvec;
    else return localUpdVec[direction];
}

////////////////////////////////////////////////////////////////////////////////

// END
