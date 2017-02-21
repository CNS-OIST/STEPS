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

// STEPS headers.
#include "steps/common.h"

#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/compdef.hpp"
#include "steps/tetexact/diff.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "third_party/easylogging++.h"

#include <iostream>

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

stex::Diff::Diff(ssolver::Diffdef * ddef, stex::Tet * tet)
: KProc()
, pDiffdef(ddef)
, pTet(tet)
, pUpdVec()
, pScaledDcst(0.0)
, pDcst(0.0)
, pCDFSelector()
, pNeighbCompLidx()
{
    assert(pDiffdef != 0);
    assert(pTet != 0);
    stex::Tet * next[4] =
    {
        pTet->nextTet(0),
        pTet->nextTet(1),
        pTet->nextTet(2),
        pTet->nextTet(3)
    };

    ligGIdx = pDiffdef->lig();
    ssolver::Compdef * cdef = pTet->compdef();
    lidxTet = cdef->specG2L(ligGIdx);

    for (uint i = 0; i < 4; ++i)
    {
        pDiffBndDirection[i] = pTet->getDiffBndDirection(i);
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
                if (pDiffBndActive[i])
                    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == cdef)
                    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0.0;
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
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
        pCDFSelector[2] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
        pCDFSelector[2] = pCDFSelector[1] + (d[2] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::Diff::~Diff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::checkpoint(std::fstream & cp_file)
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
    cp_file.write((char*)pCDFSelector, sizeof(double) * 3);
    cp_file.write((char*)pDiffBndActive, sizeof(bool) * 4);
    cp_file.write((char*)pDiffBndDirection, sizeof(bool) * 4);
    cp_file.write((char*)pNeighbCompLidx, sizeof(int) * 4);

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::restore(std::fstream & cp_file)
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
    cp_file.read((char*)pCDFSelector, sizeof(double) * 3);
    cp_file.read((char*)pDiffBndActive, sizeof(bool) * 4);
    cp_file.read((char*)pDiffBndDirection, sizeof(bool) * 4);
    cp_file.read((char*)pNeighbCompLidx, sizeof(int) * 4);

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setupDeps(void)
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
    std::set<stex::KProc*> local;

    KProcPVecCI kprocend = pTet->kprocEnd();
    for (KProcPVecCI k = pTet->kprocBegin(); k != kprocend; ++k)
    {
        // Check locally.
        if ((*k)->depSpecTet(ligGIdx, pTet) == true) {
            local.insert(*k);
        }
    }
    // Check the neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        stex::Tri * next = pTet->nextTri(i);
        if (next == 0) continue;
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTet(ligGIdx, pTet) == true) {
                local.insert(*k);
            }
        }
    }

    // Search for dependencies in neighbouring tetrahedra.
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        stex::Tet * next = pTet->nextTet(i);
        if (next == 0)
            continue;
        if (pTet->nextTri(i) != 0)
            continue;

        // Copy local dependencies.
        std::set<stex::KProc*> local2(local.begin(), local.end());

        // Find the ones 'locally' in the next tet.
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTet(ligGIdx, next) == true) {
                local2.insert(*k);
            }
        }

        // Find deps in neighbouring triangles in the next tet.
        // As said before, this cannot logically include the shared
        // triangle.
        for (uint j = 0; j < 4; ++j)
        {
            // Fetch next triangle, if it exists.
            stex::Tri * next2 = next->nextTri(j);
            if (next2 == 0) continue;

            // Find deps.
            kprocend = next2->kprocEnd();
            for (KProcPVecCI k = next2->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecTet(ligGIdx, next) == true) {
                    local2.insert(*k);
                }
            }
        }

        // Copy the set to the update vector.
        pUpdVec[i].assign(local2.begin(), local2.end());
    }
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Diff::depSpecTet(uint gidx, stex::WmVol * tet)
{
    if (pTet != tet) return false;
    if (gidx != ligGIdx) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::Diff::depSpecTri(uint gidx, stex::Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::reset(void)
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
    crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDiffBndActive(uint i, bool active)
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

bool stex::Diff::getDiffBndActive(uint i) const
{
    assert (i < 4);
    assert(pDiffBndDirection[i] == true);

    return pDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

double stex::Diff::dcst(int direction)
{
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    }
    else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDcst(double dcst)
{
    assert(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();
    
    stex::Tet * next[4] =
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
                if (pDiffBndActive[i])
                    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == pTet->compdef())
                    d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
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
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
        pCDFSelector[2] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
        pCDFSelector[2] = pCDFSelector[1] + (d[2] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Diff::setDirectionDcst(int direction, double dcst)
{
    assert(direction < 4);
    assert(direction >= 0);
    assert(dcst >= 0.0);
    directionalDcsts[direction] = dcst;
    
    stex::Tet * next[4] =
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
                if (pDiffBndActive[i])
                {
                    if (directionalDcsts.find(i) != directionalDcsts.end())
                        d[i] = (pTet->area(i) * directionalDcsts[i]) / (pTet->vol() * dist);
                    else
                        d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
                }
                else d[i] = 0.0;
            }
            else
            {
                // Bugfix 5/4/2012 IH: neighbouring tets in different compartments
                // NOT separated by a patch were allowing diffusion, even without diffusion boundary
                if (next[i]->compdef() == pTet->compdef())
                {
                    if (directionalDcsts.find(i) != directionalDcsts.end())
                        d[i] = (pTet->area(i) * directionalDcsts[i]) / (pTet->vol() * dist);
                    else
                        d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
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
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
        pCDFSelector[2] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
        pCDFSelector[2] = pCDFSelector[1] + (d[2] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

double stex::Diff::rate(steps::tetexact::Tetexact * solver)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTet->pools()[lidxTet]);
    assert(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::Diff::apply(steps::rng::RNG * rng, double dt, double simtime)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    uint * local = pTet->pools() + lidxTet;
    bool clamped = pTet->clamped(lidxTet);

    if (clamped == false)
    {
        assert(*local > 0);
    }

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    int iSel = 0;
    for (; iSel < 3; ++iSel)
        if(sel < pCDFSelector[iSel])
            break;

    // Direction iSel.
    stex::Tet * nexttet = pTet->nextTet(iSel);
    // If there is no next tet 0, pCDFSelector[0] should be zero
    // So we can assert that nextet 0 does indeed exist
    assert (nexttet != 0);
    assert(pNeighbCompLidx[iSel] > -1);

    if (nexttet->clamped(pNeighbCompLidx[iSel]) == false)
        nexttet->incCount(pNeighbCompLidx[iSel],1);

    if (clamped == false)
        pTet->incCount(lidxTet, -1);

    rExtent++;

    return pUpdVec[iSel];
}

////////////////////////////////////////////////////////////////////////////////

uint stex::Diff::updVecSize(void) const
{
    uint maxsize = pUpdVec[0].size();
    for (uint i=1; i <= 3; ++i)
    {
        if (pUpdVec[i].size() > maxsize)
            maxsize = pUpdVec[i].size();
    }
    return maxsize;
}

////////////////////////////////////////////////////////////////////////////////

// END
