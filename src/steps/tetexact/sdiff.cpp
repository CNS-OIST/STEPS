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
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"
#include "steps/solver/diffdef.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/tetexact/sdiff.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/tetexact.hpp"

// third party headers
#include "third_party/easylogging++.h"



////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace smath = steps::math;

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::SDiff(ssolver::Diffdef * sdef, stex::Tri * tri)
: KProc()
, pSDiffdef(sdef)
, pTri(tri)
, pUpdVec()
, pScaledDcst(0.0)
, pDcst(0.0)
, pCDFSelector()

{
    assert(pSDiffdef != 0);
    assert(pTri != 0);
    stex::Tri * next[3] =
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
        pDiffBndDirection[i] = pTri->getDiffBndDirection(i);
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
            if (pDiffBndDirection[i] == true)
            {
                if (pDiffBndActive[i])
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
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

stex::SDiff::~SDiff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::checkpoint(std::fstream & cp_file)
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
    cp_file.write((char*)pCDFSelector, sizeof(double) * 2);
    cp_file.write((char*)pDiffBndActive, sizeof(bool) * 3);
    cp_file.write((char*)pDiffBndDirection, sizeof(bool) * 3);
    cp_file.write((char*)pNeighbPatchLidx, sizeof(int) * 3);

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::restore(std::fstream & cp_file)
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
    cp_file.read((char*)pCDFSelector, sizeof(double) * 2);
    cp_file.read((char*)pDiffBndActive, sizeof(bool) * 3);
    cp_file.read((char*)pDiffBndDirection, sizeof(bool) * 3);
    cp_file.read((char*)pNeighbPatchLidx, sizeof(int) * 3);

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setupDeps(void)
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


    // Search for dependencies in the 'source' triangle.
    std::set<stex::KProc*> local;

    KProcPVecCI kprocend = pTri->kprocEnd();
    for (KProcPVecCI k = pTri->kprocBegin(); k != kprocend; ++k)
    {
        // Check locally.
        if ((*k)->depSpecTri(ligGIdx, pTri) == true) {
            local.insert(*k);
        }
    }

    // Check the neighbouring tetrahedrons.
    stex::WmVol * itet[2] = {pTri->iTet(), pTri->oTet()};
    for (uint i = 0; i < 2; ++i)
    {
        if (itet[i] == 0)
            continue;
        kprocend = itet[i]->kprocEnd();
        for (KProcPVecCI k = itet[i]->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTri(ligGIdx, pTri) == true) {
                local.insert(*k);
            }
        }
    }

    // Search for dependencies in neighbouring triangles.
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangle, if it exists.
        stex::Tri * next = pTri->nextTri(i);
        if (next == 0) continue;

        // Copy local dependencies.
        std::set<stex::KProc*> local2(local.begin(), local.end());

        // Find the ones 'locally' in the next tri.
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTri(ligGIdx, next) == true) {
                local2.insert(*k);
            }
        }

        // Fetch inner tetrahedron, if it exists.
        stex::WmVol * itet[2] = {pTri->iTet(), pTri->oTet()};
        for (uint j = 0; j < 2; ++j)
        {
            if (itet[j] == 0)
                continue;

            // Find deps.
            kprocend = itet[j]->kprocEnd();
            for (KProcPVecCI k = itet[j]->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecTri(ligGIdx, next) == true) {
                    local2.insert(*k);
                }
            }
        }

        // Copy the set to the update vector.
        pUpdVec[i].assign(local2.begin(), local2.end());
    }

}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTet(uint gidx, stex::WmVol * tet)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::SDiff::depSpecTri(uint gidx, stex::Tri * tri)
{
    if (pTri != tri) return false;
    if (gidx != ligGIdx) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::reset(void)
{
    resetExtent();

    pDiffBndActive[0] = false;
    pDiffBndActive[1] = false;
    pDiffBndActive[2] = false;

    uint ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;

}

////////////////////////////////////////////////////////////////////////////////

double stex::SDiff::dcst(int direction)
{
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    }
    else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setDcst(double dcst)
{
    assert(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    stex::Tri * next[3] =
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
            if (pDiffBndDirection[i] == true)
            {
                if (pDiffBndActive[i])
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
    assert(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::SDiff::setDirectionDcst(int direction, double dcst)
{
    #ifdef DIRECTIONAL_DCST_DEBUG
    CLOG(DEBUG, "steps_debug") << "Direction: " << direction << " dcst: " << dcst << "\n";
    #endif
    
    assert(direction < 3);
    assert(direction >= 0);
    assert(dcst >= 0.0);
    directionalDcsts[direction] = dcst;
    
    stex::Tri * next[3] =
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
            if (pDiffBndDirection[i] == true)
            {
                if (pDiffBndActive[i])
                {
                    if (directionalDcsts.find(i) != directionalDcsts.end())
                        d[i] = (pTri->length(i) * directionalDcsts[i]) / (pTri->area() * dist);
                    else
                        d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
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
                        d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
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
    assert(pScaledDcst >= 0);
    
    // Setup the selector distribution.
    if (pScaledDcst == 0.0)
    {
        pCDFSelector[0] = 0.0;
        pCDFSelector[1] = 0.0;
    }
    else
    {
        pCDFSelector[0] = d[0] / pScaledDcst;
        pCDFSelector[1] = pCDFSelector[0] + (d[1] / pScaledDcst);
    }
}
////////////////////////////////////////////////////////////////////////////////

double stex::SDiff::rate(steps::tetexact::Tetexact * solver)
{
    if (inactive()) return 0.0;

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTri->pools()[lidxTri]);
    assert(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::SDiff::apply(steps::rng::RNG * rng, double dt, double simtime)
{
    //uint lidxTet = this->lidxTet;
    // Pre-fetch some general info.


    // Apply local change.
    uint * local = pTri->pools() + lidxTri;
    bool clamped = pTri->clamped(lidxTri);

    if (clamped == false)
    {
        assert(*local > 0);
    }

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    int iSel = 0;
    for (; iSel < 2; ++iSel)
        if(sel < pCDFSelector[iSel])
            break;

    // Direction iSel.
    stex::Tri * nexttri = pTri->nextTri(iSel);
    // If there is no next tet iSel, pCDFSelector[iSel] should be zero
    // So we can assert that nextet 0 does indeed exist
    assert (nexttri != 0);

    if (nexttri->clamped(lidxTri) == false)
        nexttri->incCount(lidxTri,1);

    if (clamped == false)
        pTri->incCount(lidxTri, -1);

    rExtent++;

    return pUpdVec[iSel];
}

////////////////////////////////////////////////////////////////////////////////

uint stex::SDiff::updVecSize(void) const
{
    uint maxsize = pUpdVec[0].size();
    for (uint i=1; i <= 2; ++i)
    {
        if (pUpdVec[i].size() > maxsize)
            maxsize = pUpdVec[i].size();
    }
    return maxsize;
}

////////////////////////////////////////////////////////////////////////////////

// END
