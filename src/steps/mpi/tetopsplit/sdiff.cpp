/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
#include <iostream>
#include <fstream>
#include <sstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/mpi/tetopsplit/sdiff.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

smtos::SDiff::SDiff(ssolver::Diffdef * sdef, smtos::Tri * tri)
: KProc()
, pSDiffdef(sdef)
, pTri(tri)
, localAllUpdVec()
, remoteAllUpdVec()
, pEmptyvec()
, idxEmptyvec()
, pScaledDcst(0.0)
, pDcst(0.0)
, pNonCDFSelector()
, pDirections()
, pNdirections(0)
, pSDiffBndActive()
, pSDiffBndDirection()
{
    AssertLog(pSDiffdef != 0);
    AssertLog(pTri != 0);
    type = KP_SDIFF;
    smtos::Tri * next[3] =
    {
        pTri->nextTri(0),
        pTri->nextTri(1),
        pTri->nextTri(2),
    };
    
    ligGIdx = pSDiffdef->lig();
    ssolver::Patchdef * pdef = pTri->patchdef();
    lidxTri = pdef->specG2L(ligGIdx);

    for (uint i = 0; i < 3; ++i)
    {
        pSDiffBndDirection[i] = pTri->getSDiffBndDirection(i);
        if (next[i] == 0)
        {
            pNeighbPatchLidx[i] = -1;
            continue;
        }
        else
        {
            pNeighbPatchLidx[i] = next[i]->patchdef()->specG2L(ligGIdx);
        }
    }

    // Precalculate part of the scaled diffusion constant.
    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    pDcst = dcst;

    double d[3] = { 0.0, 0.0, 0.0};
    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pSDiffBndDirection[i] == true)
            {
                if (pSDiffBndActive[i])
                    d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                else d[i] = 0.0;
            }
            else
            {
                if (next[i]->patchdef() == pTri->patchdef())
                    d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                else d[i] = 0.0;
            }
        }
    }
    // Compute scaled "diffusion constant".
    for (uint i = 0; i < 3; ++i)
    {
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;

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
    }
}

////////////////////////////////////////////////////////////////////////////////

smtos::SDiff::~SDiff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));

    uint n_direct_dcsts = directionalDcsts.size();
    cp_file.write((char*)&n_direct_dcsts, sizeof(uint));
    for (auto& item : directionalDcsts) {
        cp_file.write((char*)&(item.first), sizeof(uint));
        cp_file.write((char*)&(item.second), sizeof(double));
    }
    
    cp_file.write((char*)&pScaledDcst, sizeof(double));
    cp_file.write((char*)&pDcst, sizeof(double));

    cp_file.write((char*)pNonCDFSelector, sizeof(double) * 3);
    cp_file.write((char*)pSDiffBndActive, sizeof(bool) * 3);
    cp_file.write((char*)pSDiffBndDirection, sizeof(bool) * 3);
    cp_file.write((char*)pNeighbPatchLidx, sizeof(int) * 3);

    // Need to add directional stuff here if checkpointing ever implemented

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::restore(std::fstream & cp_file)
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
    
    cp_file.read((char*)&pScaledDcst, sizeof(double));
    cp_file.read((char*)&pDcst, sizeof(double));

    cp_file.read((char*)pNonCDFSelector, sizeof(double) * 3);
    cp_file.read((char*)pSDiffBndActive, sizeof(bool) * 3);
    cp_file.read((char*)pSDiffBndDirection, sizeof(bool) * 3);
    cp_file.read((char*)pNeighbPatchLidx, sizeof(int) * 3);

    // Need to add directional stuff here if checkpointing ever implemented
    
    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::setupDeps(void)
{
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' triangle
    //   * any neighbouring tetrahedrons
    //
    // But also in the possible 'destination' triangles (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring triangles
    //   * any neighbouring tetrahedrons of these neighbouring tris
    //
    AssertLog(pTri->getInHost());

    // Search for dependencies in the 'source' triangle.
    std::set<uint> remote;
    std::set<uint> remote_all;
    
    std::set<smtos::KProc*> local;
    std::set<smtos::KProc*> local_all;
    
    uint nkprocs = pTri->countKProcs();
    uint startKProcIdx = pTri->getStartKProcIdx();
    for (uint sk = 0; sk < nkprocs; sk++)
    {
        // Check locally.
        if (pTri->KProcDepSpecTri(sk, pTri, ligGIdx) == true) {
            local.insert(pTri->getKProc(sk));
        }
    }

    // Check the neighbouring tetrahedrons.
    smtos::WmVol * itet = pTri->iTet();
    if (itet != 0)
    {
        if (pTri->getHost() != itet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }
        nkprocs = itet->countKProcs();
        startKProcIdx = itet->getStartKProcIdx();
        for (uint k = 0; k < nkprocs; k++)
        {
            if (itet->KProcDepSpecTri(k, pTri, ligGIdx) == true) {
                local.insert(itet->getKProc(k));
            }
        }
    }

    smtos::WmVol * otet = pTri->oTet();
    if (otet != 0)
    {
        if (pTri->getHost() != otet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << otet->idx()  << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }
        nkprocs = otet->countKProcs();
        startKProcIdx = otet->getStartKProcIdx();
        for (uint k = 0; k < nkprocs; k++)
        {
            if (otet->KProcDepSpecTri(k, pTri, ligGIdx) == true) {
                local.insert(otet->getKProc(k));
            }
        }
    }


    // Search for dependencies in neighbouring triangles.
    // which can be in different host
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangle, if it exists.
        smtos::Tri * next = pTri->nextTri(i);
        if (next == 0) continue;
        
        if (next->getHost() != pTri->getHost()) {
            pTri->solver()->addNeighHost(next->getHost());
            pTri->solver()->registerBoundaryTri(next);
        }
        
        // Copy local dependencies.
        std::set<smtos::KProc*> local2(local.begin(), local.end());
        std::set<uint> remote2;

        // Find the ones 'locally' in the next tri.
        nkprocs = next->countKProcs();
        startKProcIdx = next->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            // Check locally.
            if (next->KProcDepSpecTri(sk, next, ligGIdx) == true) {
                if (next->getHost() == pTri->getHost()) {
                    local2.insert(next->getKProc(sk));
                }
                else {
                    remote2.insert(startKProcIdx + sk);
                }
            }
        }
        // Fetch inner tetrahedron, if it exists.
        smtos::WmVol * itet = next->iTet();
        if (itet != 0)
        {
            if (next->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            
            // Find deps.
            nkprocs = itet->countKProcs();
            startKProcIdx = itet->getStartKProcIdx();
            for (uint k = 0; k < nkprocs; k++)
            {
                if (itet->KProcDepSpecTri(k, next, ligGIdx) == true) {
                    if (itet->getHost() == pTri->getHost()) {
                        local2.insert(itet->getKProc(k));
                    }
                    else {
                        remote2.insert(startKProcIdx + k);
                    }
                }
            }
        }

        // Fetch outer tetrahedron, if it exists.
        smtos::WmVol * otet = next->oTet();
        if (otet != 0)
        {
            if (next->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next->idx() << " and its compartment tetrahedron " << otet->idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            // Find deps.
            nkprocs = otet->countKProcs();
            startKProcIdx = otet->getStartKProcIdx();
            for (uint k = 0; k < nkprocs; k++)
            {
                if (otet->KProcDepSpecTri(k, next, ligGIdx) == true) {
                    if (otet->getHost() == pTri->getHost()) {
                        local2.insert(otet->getKProc(k));
                    }
                    else {
                        remote2.insert(startKProcIdx + k);
                    }
                }
            }
        }

        localUpdVec[i].assign(local2.begin(), local2.end());
        remoteUpdVec[i].assign(remote2.begin(), remote2.end());

        local_all.insert(local2.begin(), local2.end());
        remote_all.insert(remote2.begin(), remote2.end());

    }
    localAllUpdVec.assign(local_all.begin(), local_all.end());
    remoteAllUpdVec.assign(remote_all.begin(), remote_all.end());
    
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::SDiff::depSpecTet(uint gidx, smtos::WmVol * tet)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::SDiff::depSpecTri(uint gidx, smtos::Tri * tri)
{
    if (pTri != tri) return false;
    if (gidx != ligGIdx) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::reset(void)
{
    resetExtent();

    pSDiffBndActive[0] = false;
    pSDiffBndActive[1] = false;
    pSDiffBndActive[2] = false;

    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    // cannot reset because the position is used in solver
    //crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

double smtos::SDiff::dcst(int direction)
{
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    }
    else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::setDcst(double dcst)
{
    AssertLog(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    smtos::Tri * next[3] =
    {
        pTri->nextTri(0),
        pTri->nextTri(1),
        pTri->nextTri(2),
    };

    double d[3] = { 0.0, 0.0, 0.0};

    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pSDiffBndDirection[i] == true)
            {
                if (pSDiffBndActive[i])
                    d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                else d[i] = 0;
            }
            else
            {
                if (next[i]->patchdef() == pTri->patchdef())
                    d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                else d[i] = 0;
            }
        }
    }

    // Compute scaled "diffusion constant".
    pScaledDcst = 0.0;

    for (uint i = 0; i < 3; ++i)
    {
        pScaledDcst += d[i];
    }
    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    pNdirections = 0;
    pDirections.clear();
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;

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
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::setDirectionDcst(int direction, double dcst)
{
    AssertLog(direction < 3);
    AssertLog(direction >= 0);
    AssertLog(dcst >= 0.0);
    directionalDcsts[direction] = dcst;
    
    // Automatically activate boundary diffusion if necessary
    if (pSDiffBndDirection[direction] == true) pSDiffBndActive[direction] = true;

    smtos::Tri * next[3] =
    {
        pTri->nextTri(0),
        pTri->nextTri(1),
        pTri->nextTri(2)
    };
    
    double d[3] = { 0.0, 0.0, 0.0};

    
    for (uint i = 0; i < 3; ++i)
    {
        // if directional diffusion dcst exists use directional dcst, else use the standard one
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != 0))
        {
            if (pSDiffBndDirection[i] == true)
            {
                if (pSDiffBndActive[i])
                {
                    if (directionalDcsts.find(i) != directionalDcsts.end())
                        d[i] = (pTri->length(i) * directionalDcsts[i]) / (pTri->area() * dist);
                    else
                    	// This part must use the default pDcst, not the function argument dcst
                        d[i] = (pTri->length(i) * pDcst) / (pTri->area() * dist);
                }
                else d[i] = 0;
            }
            else
            {
                if (next[i]->patchdef() == pTri->patchdef())
                {
                    if (directionalDcsts.find(i) != directionalDcsts.end())
                        d[i] = (pTri->length(i) * directionalDcsts[i]) / (pTri->area() * dist);
                    else
                    	// This part must use the default pDcst, not the function argument dcst
                        d[i] = (pTri->length(i) * pDcst) / (pTri->area() * dist);
                }
                else d[i] = 0;
            }
        }
    }
    
    // Compute scaled "diffusion constant".
    pScaledDcst = 0.0;

    for (uint i = 0; i < 3; ++i)
    {
        pScaledDcst += d[i];
    }
    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    pNdirections = 0;
    pDirections.clear();
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector[0] = 0.0;
        pNonCDFSelector[1] = 0.0;
        pNonCDFSelector[2] = 0.0;
    }
    else
    {
        pNonCDFSelector[0] = d[0]/pScaledDcst;
        pNonCDFSelector[1] = d[1]/pScaledDcst;
        pNonCDFSelector[2] = d[2]/pScaledDcst;

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
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::SDiff::rate(steps::mpi::tetopsplit::TetOpSplitP * solver)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTri->pools()[lidxTri]);
    AssertLog(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::SDiff::getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver)
{
    return pScaledDcst;
}

////////////////////////////////////////////////////////////////////////////////

int smtos::SDiff::apply(steps::rng::RNG * rng)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    uint * local = pTri->pools() + lidxTri;
    bool clamped = pTri->clamped(lidxTri);

	if (clamped == false)
	{
		//AssertLog(*local > 0);
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
        smtos::Tri * nexttri = pTri->nextTri(0);
        // If there is no next tet 0, pCDFSelector[0] should be zero
        // So we can assert that nextet 0 does indeed exist
        AssertLog(nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
            nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return 0;
    }
    else if (sel < pNonCDFSelector[0] + pNonCDFSelector[1])
    {
        // Direction 2.
        smtos::Tri * nexttri = pTri->nextTri(1);
        // If there is no next tet 1, pCDFSelector[1] should be zero
        // So we can assert that nextet 1 does indeed exist
        AssertLog(nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
            nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return 1;
    }
    else
    {
        // Direction 3.
        smtos::Tri * nexttri = pTri->nextTri(2);
        AssertLog(nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
        {
            nexttri->incCount(lidxTri,1);
        }
        if (clamped == false) {pTri->incCount(lidxTri, -1); }

        rExtent++;

        return 2;

    }

    // This should never happen!
    ProgErrLog("Cannot find a suitable direction for surface diffusion!");
}

///////////////////////////////////////////////////////////////////////////////

int smtos::SDiff::apply(steps::rng::RNG * rng, uint nmolcs)
{
	// Apply local change.
    uint * local = pTri->pools() + lidxTri;
	bool clamped = pTri->clamped(lidxTri);

	if (clamped == false)
	{
		//AssertLog(*local > 0);
		if (*local == 0)
		{
			return -2;
		}
	}

	AssertLog(pNdirections >= 1);

    // Multinomial
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

        	smtos::Tri * nexttri = pTri->nextTri(direction);

            AssertLog(nexttri != 0);

            if (nexttri->clamped(lidxTri) == false)
        	{
            	nexttri->incCount(lidxTri,molcsthisdir);
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
    	smtos::Tri * nexttri = pTri->nextTri(direction);

        AssertLog(nexttri != 0);

        if (nexttri->clamped(lidxTri) == false)
    	{
        	nexttri->incCount(lidxTri,molcsthisdir);
    	}

    	molcs_moved+=molcsthisdir;
    }


	AssertLog(molcs_moved == nmolcs);


	if (clamped == false) {pTri->incCount(lidxTri, -nmolcs); }

	rExtent+=nmolcs;

	return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::SDiff::getRemoteUpdVec(int direction)
{
    if (direction == -1) return remoteAllUpdVec;
    else if (direction == -2) return idxEmptyvec;
    else return remoteUpdVec[direction];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::SDiff::getLocalUpdVec(int direction)
{
    if (direction == -1) return localAllUpdVec;
    else if (direction == -2) return pEmptyvec;
    else return localUpdVec[direction];
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::setSDiffBndActive(uint i, bool active)
{
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    // Only need to update if the flags are changing
    if (pSDiffBndActive[i] != active)
    {
        pSDiffBndActive[i] = active;
        setDcst(pDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::SDiff::getSDiffBndActive(uint i) const
{
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    return pSDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

// END
