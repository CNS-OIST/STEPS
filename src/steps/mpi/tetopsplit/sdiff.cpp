/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/sdiff.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/patchdef.hpp"

// logging
#include "easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

smtos::SDiff::SDiff(ssolver::Diffdef * sdef, smtos::Tri * tri)
: 
 pSDiffdef(sdef)
, pTri(tri)
{
    AssertLog(pSDiffdef != nullptr);
    AssertLog(pTri != nullptr);
    type = KP_SDIFF;
    std::array<smtos::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };

    ssolver::Patchdef * pdef = pTri->patchdef();
    lidxTri = pdef->specG2L(pSDiffdef->lig());

    // Precalculate part of the scaled diffusion constant.
    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    pDcst = dcst;

    std::array<double, 3> d{0.0, 0.0, 0.0};

    for (uint i = 0; i < 3; ++i)
    {
        pSDiffBndDirection[i] = pTri->getSDiffBndDirection(i);
        if (next[i] != nullptr) {
            pNeighbPatchLidx[i] = next[i]->patchdef()->specG2L(pSDiffdef->lig());
        }

        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef()) {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
                pScaledDcst += d[i];
            }
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst > 0.0)
    {
        pNonCDFSelector = { d[0]/pScaledDcst, 
                            d[1]/pScaledDcst, 
                            d[2]/pScaledDcst
                            };

        for (uint i = 0; i < 3; ++i) { 
            if (d[i] > 0.0)
            {
			    pDirections.push_back(i);
			    pNdirections+=1;
       	    }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    auto n_direct_dcsts = directionalDcsts.size();
    cp_file.write(reinterpret_cast<char*>(&n_direct_dcsts), sizeof(uint));
    for (auto& item : directionalDcsts) {
        cp_file.write(reinterpret_cast<const char*>(&item.first), sizeof(uint));
        cp_file.write(reinterpret_cast<char*>(&item.second), sizeof(double));
    }
    
    cp_file.write(reinterpret_cast<char*>(&pScaledDcst), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pDcst), sizeof(double));

    cp_file.write(reinterpret_cast<char*>(pNonCDFSelector.data()), sizeof(double) * 3);
    cp_file.write(reinterpret_cast<char*>(pSDiffBndActive.data()), sizeof(bool) * 3);
    cp_file.write(reinterpret_cast<char*>(pSDiffBndDirection.data()), sizeof(bool) * 3);
    cp_file.write(reinterpret_cast<char*>(pNeighbPatchLidx.data()), sizeof(ssolver::lidxT) * 3);

    // Need to add directional stuff here if checkpointing ever implemented

    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.read(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    uint n_direct_dcsts = 0;
    cp_file.read(reinterpret_cast<char*>(&n_direct_dcsts), sizeof(uint));
    for (uint i = 0; i < n_direct_dcsts; i++) {
        uint id = 0;
        double value = 0.0;
        cp_file.read(reinterpret_cast<char*>(&id), sizeof(uint));
        cp_file.read(reinterpret_cast<char*>(&value), sizeof(double));
        directionalDcsts[id] = value;
    }
    
    cp_file.read(reinterpret_cast<char*>(&pScaledDcst), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pDcst), sizeof(double));

    cp_file.read(reinterpret_cast<char*>(pNonCDFSelector.data()), sizeof(double) * 3);
    cp_file.read(reinterpret_cast<char*>(pSDiffBndActive.data()), sizeof(bool) * 3);
    cp_file.read(reinterpret_cast<char*>(pSDiffBndDirection.data()), sizeof(bool) * 3);
    cp_file.read(reinterpret_cast<char*>(pNeighbPatchLidx.data()), sizeof(ssolver::lidxT) * 3);

    // Need to add directional stuff here if checkpointing ever implemented
    
    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::setupDeps()
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
    for (uint sk = 0; sk < nkprocs; sk++)
    {
        // Check locally.
        if (pTri->KProcDepSpecTri(sk, pTri, pSDiffdef->lig())) {
            local.insert(pTri->getKProc(sk));
        }
    }

    // Check the neighbouring tetrahedrons.
    {
        smtos::WmVol * itet = pTri->iTet();
        if (itet != nullptr)
        {
            if (pTri->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = itet->countKProcs();
            for (uint k = 0; k < nkprocs; k++)
            {
                if (itet->KProcDepSpecTri(k, pTri, pSDiffdef->lig())) {
                    local.insert(itet->getKProc(k));
                }
            }
        }
    }

    {
        smtos::WmVol *otet = pTri->oTet();
        if (otet != nullptr) {
            if (pTri->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << otet->idx()
                   << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = otet->countKProcs();
            for (uint k = 0; k < nkprocs; k++) {
                if (otet->KProcDepSpecTri(k, pTri, pSDiffdef->lig())) {
                    local.insert(otet->getKProc(k));
                }
            }
        }
    }


    // Search for dependencies in neighbouring triangles.
    // which can be in different host
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangle, if it exists.
        smtos::Tri * next = pTri->nextTri(i);
        if (next == nullptr) { continue;
}
        
        if (next->getHost() != pTri->getHost()) {
            pTri->solver()->addNeighHost(next->getHost());
            pTri->solver()->registerBoundaryTri(next);
        }
        
        // Copy local dependencies.
        std::set<smtos::KProc*> local2(local.begin(), local.end());
        std::set<uint> remote2;

        // Find the ones 'locally' in the next tri.
        nkprocs = next->countKProcs();
        auto startKProcIdx = next->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            // Check locally.
            if (next->KProcDepSpecTri(sk, next, pSDiffdef->lig())) {
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
        if (itet != nullptr)
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
                if (itet->KProcDepSpecTri(k, next, pSDiffdef->lig())) {
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
        if (otet != nullptr)
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
                if (otet->KProcDepSpecTri(k, next, pSDiffdef->lig())) {
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

bool smtos::SDiff::depSpecTet(uint /*gidx*/, smtos::WmVol * /*tet*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::SDiff::depSpecTri(uint gidx, smtos::Tri * tri)
{
    if (pTri != tri) { return false;
}
    if (gidx != pSDiffdef->lig()) { return false;
}
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::SDiff::reset()
{
    resetExtent();

    pSDiffBndActive = {false, false, false};

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
    auto search_result = directionalDcsts.find(direction);
    if (search_result != directionalDcsts.end()) {
        return search_result->second;
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

    std::array<smtos::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };

    std::array<double, 3> d{0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    for (uint i = 0; i < 3; ++i)
    {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) || (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef()))
            {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
            }
        }
        pScaledDcst += d[i];
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    pNdirections = 0;
    pDirections.clear();
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector = {0.0, 0.0, 0.0};
    }
    else
    {
        pNonCDFSelector = 
        { 
            d[0]/pScaledDcst, 
            d[1]/pScaledDcst, 
            d[2]/pScaledDcst
        };

        for (uint i = 0; i < 3; ++i) { 
            if (d[i] > 0.0)
            {
			    pDirections.push_back(i);
			    pNdirections+=1;
       	    }
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
    if (pSDiffBndDirection[direction] == true) { 
        pSDiffBndActive[direction] = true;
    }

    std::array<smtos::Tri *, 3> next{pTri->nextTri(0),
                                    pTri->nextTri(1),
                                    pTri->nextTri(2)
                                    };
    
    std::array<double, 3> d{0.0, 0.0, 0.0};
    pScaledDcst = 0.0;
    
    for (uint i = 0; i < 3; ++i)
    {
        // if directional diffusion dcst exists use directional dcst, else use the standard one
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr))
        {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) || (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef())) {
                auto search_result = directionalDcsts.find(i);
                if (search_result != directionalDcsts.end()) {
                    d[i] = (pTri->length(i) * search_result->second) / (pTri->area() * dist);
                }
                else {
                    // This part must use the default pDcst, not the function argument dcst
                    d[i] = (pTri->length(i) * pDcst) / (pTri->area() * dist);
                }
            }
        }
        pScaledDcst += d[i];
    }
    
    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    pNdirections = 0;
    pDirections.clear();
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pNonCDFSelector = {0.0, 0.0, 0.0};
    }
    else
    {
        pNonCDFSelector = 
        { 
            d[0]/pScaledDcst, 
            d[1]/pScaledDcst, 
            d[2]/pScaledDcst
        };

        for (uint i = 0; i < 3; ++i) { 
            if (d[i] > 0.0)
            {
			    pDirections.push_back(i);
			    pNdirections+=1;
       	    }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double smtos::SDiff::rate(steps::mpi::tetopsplit::TetOpSplitP * /*solver*/)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTri->pools()[lidxTri]);
    AssertLog(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::SDiff::getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * /*solver*/) const
{
    return pScaledDcst;
}

////////////////////////////////////////////////////////////////////////////////

int smtos::SDiff::apply(const rng::RNGptr &rng)
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

    int iSel = 0;
    double pCDFSelector = 0.0;
    for (; iSel < 2; ++iSel) {
        pCDFSelector += pNonCDFSelector[iSel];
        if(sel < pCDFSelector) {
            break;
        }
    }       

    // Direction iSel.
    smtos::Tri * nexttri = pTri->nextTri(iSel);
    AssertLog(nexttri != nullptr);

    if (nexttri->clamped(pNeighbPatchLidx[iSel]) == false) {
        nexttri->incCount(pNeighbPatchLidx[iSel],1);
}

    if (clamped == false) {
        pTri->incCount(lidxTri, -1);
}

    rExtent++;
    return iSel;
}

///////////////////////////////////////////////////////////////////////////////

int smtos::SDiff::apply(const rng::RNGptr &rng, uint nmolcs)
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
        if (chance >= 1.0) { chance=1.0;
}

        uint molcsthisdir = rng->getBinom(max_molcs, chance);

        if (molcsthisdir != 0u)
        {

        	smtos::Tri * nexttri = pTri->nextTri(direction);

            AssertLog(nexttri != nullptr);
            AssertLog(pNeighbPatchLidx[direction] != ssolver::LIDX_UNDEFINED);

            if (nexttri->clamped(pNeighbPatchLidx[direction]) == false)
        	{
            	nexttri->incCount(pNeighbPatchLidx[direction],molcsthisdir);
        	}

        	molcs_moved+=molcsthisdir;
        }
        if (molcs_moved == nmolcs) { break;
}

    }
    //last direction, chance =1
    uint direction = pDirections[pNdirections-1];
    int molcsthisdir = nmolcs-molcs_moved;
    if (molcsthisdir != 0)
    {
    	smtos::Tri * nexttri = pTri->nextTri(direction);

        AssertLog(nexttri != nullptr);
        AssertLog(pNeighbPatchLidx[direction] != ssolver::LIDX_UNDEFINED);

        if (nexttri->clamped(pNeighbPatchLidx[direction]) == false)
    	{
        	nexttri->incCount(pNeighbPatchLidx[direction],molcsthisdir);
    	}

    	molcs_moved+=molcsthisdir;
    }


	AssertLog(molcs_moved == nmolcs);


	if (clamped == false) {pTri->incCount(lidxTri, -nmolcs); }

	rExtent+=nmolcs;

	return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::SDiff::getRemoteUpdVec(int direction) const
{
    if (direction == -1) return remoteAllUpdVec;
    else if (direction == -2) return idxEmptyvec;
    else return remoteUpdVec[direction];
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::SDiff::getLocalUpdVec(int direction) const
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
