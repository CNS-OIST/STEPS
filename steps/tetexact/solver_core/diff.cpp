////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
// 
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id:diff.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <algorithm>
#include <cassert>
#include <cmath>
#include <iostream>
#include <iterator>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(ssim::DiffDef * ddef, Tet * tet)
: KProc()
, pDiffDef(ddef)
, pTet(tet)
, pUpdVec()
, pScaledDcst(0.0)
, pCDFSelector()
{
    // Fetch neighbouring voxels.
    Tet * next[4] = 
    { 
        pTet->nextTet(0), 
        pTet->nextTet(1), 
        pTet->nextTet(2), 
        pTet->nextTet(3) 
    };
    
    // Precalculate part of the scaled diffusion constant.
    double dcst = pDiffDef->dcst();
    double d[4] = { 0.0, 0.0, 0.0, 0.0 };
    for (uint i = 0; i < 4; ++i)
    {
        // Compute the scaled diffusion constant.
        double dist = pTet->dist(i);
        if (dist > 0.0)
            d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
    }
    
    // Compute scaled "diffusion constant".
    pScaledDcst = d[0] + d[1] + d[2] + d[3];
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

Diff::~Diff(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setupDeps(void)
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
    
    uint gidx = def()->lig();
    
    // Search for dependencies in the 'source' tetrahedron.
    SchedIDXVec local;
    KProcPVecCI kprocend = pTet->kprocEnd();
    for (KProcPVecCI k = pTet->kprocBegin(); k != kprocend; ++k)
    {
        // Check locally.
        if ((*k)->depSpecTet(gidx, pTet) == true)
            local.push_back((*k)->schedIDX());
    }
    // Check the neighbouring triangles.
    for (uint i = 0; i < 4; ++i)
    {
        Tri * next = pTet->nextTri(i);
        if (next == 0) continue;
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTet(gidx, pTet) == true)
                pUpdVec[i].push_back((*k)->schedIDX());
        }
    }
    
    // Search for dependencies in neighbouring tetrahedrons.
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        Tet * next = pTet->nextTet(i);
        if (next == 0) continue;
        if (pTet->nextTri(i) != 0) continue;
        
        // Copy local dependencies.
        std::copy(local.begin(), local.end(), 
            std::inserter(pUpdVec[i], pUpdVec[i].end()));
        
        // Find the ones 'locally' in the next tet.
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            if ((*k)->depSpecTet(gidx, next) == true)
                pUpdVec[i].push_back((*k)->schedIDX());
        }
        
        // Find deps in neighbouring triangles in the next tet.
        // As said before, this cannot logically include the shared 
        // triangle.
        for (uint j = 0; j < 4; ++j)
        {
            // Fetch next triangle, if it exists.
            Tri * next2 = next->nextTri(j);
            if (next2 == 0) continue;
            
            // Find deps.
            kprocend = next2->kprocEnd(); 
            for (KProcPVecCI k = next2->kprocBegin(); k != kprocend; ++k)
            {
                if ((*k)->depSpecTet(gidx, next) == true)
                    pUpdVec[i].push_back((*k)->schedIDX());
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Diff::depSpecTet(uint gidx, Tet * tet)
{
    if (pTet != tet) return false;
    if (gidx != def()->lig()) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool Diff::depSpecTri(uint gidx, Tri * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::reset(void)
{
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

double Diff::rate(void) const
{
    if (inactive()) return 0.0;
    
    // Pre-fetch some general info.
    ssim::CompDef * cdef = pTet->compdef();
    // Fetch the ligand as global index.
    uint gidx = pDiffDef->lig();
    // As local index.
    uint lidx = cdef->specG2L(gidx);
    
    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTet->pools()[lidx]);
    assert(std::isnan(rate) == false);
    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

SchedIDXVec const & Diff::apply(State * s)
{
    // Pre-fetch some general info.
    ssim::CompDef * cdef = pTet->compdef();
    // Fetch the ligand as global index.
    uint gidx = def()->lig();
    // As local index.
    uint lidx = cdef->specG2L(gidx);
    
    // Apply local change.
    uint * local = pTet->pools() + lidx;
    assert(*local > 0);
    *local -= 1;
    
    // Apply change in next voxel: select a direction.
    double sel = s->rng()->getUnfEE();
    if (sel < pCDFSelector[0])
    {
        // Direction 1.
        uint * next = pTet->nextTet(0)->pools() + lidx;
        *next += 1;
        return pUpdVec[0];
    }
    else if (sel < pCDFSelector[1])
    {
        // Direction 2.
        uint * next = pTet->nextTet(1)->pools() + lidx;
        *next += 1;
        return pUpdVec[1];
    }
    else if (sel < pCDFSelector[2])
    {
        // Direction 3.
        uint * next = pTet->nextTet(2)->pools() + lidx;
        *next += 1;
        return pUpdVec[2];
    }
    else 
    {
        // Direction 4.
        uint * next = pTet->nextTet(3)->pools() + lidx;
        *next += 1;
        return pUpdVec[3];
    }
    
    // This should never happen!
    assert(0);
    return pUpdVec[0];
}

////////////////////////////////////////////////////////////////////////////////

// END
