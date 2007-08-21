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

////////////////////////////////////////////////////////////////////////////////

Diff::Diff(DiffDef * ddef, Tet * tet)
: pDiffDef(ddef)
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
	// Fetch ligand index.
	uint gidx = def()->lig();
	
	// Search for local dependencies.
	SchedIDXVec local;
	// First check for other diffusion processes.
	std::vector<Diff*>::const_iterator diffend = pTet->diffEnd();
	for (std::vector<Diff*>::const_iterator d = pTet->diffBegin(); 
		d != diffend; ++d)
	{
		if ((*d)->depSpecTet(gidx, pTet) == true)
			local.push_back((*d)->schedIDX());
	}
	
	// Search for dependencies in neighbouring tetrahedrons.
	for (uint i = 0; i < 4; ++i)
	{
		// Fetch next tetrahedron, if it exists.
		Tet * next = pTet->nextTet(i);
		if (next == 0) continue;
		// Later, also check if there is a triangle that might block
		// diffusion.
		
		// Copy local dependencies.
		std::copy(local.begin(), local.end(), 
			std::inserter(pUpdVec[i], pUpdVec[i].end()));
		
		// First check for diffusion processes.
		diffend = next->diffEnd();
		for (std::vector<Diff*>::const_iterator d = next->diffBegin(); 
			d != diffend; ++d)
		{
			if ((*d)->depSpecTet(gidx, next) == true)
				pUpdVec[i].push_back((*d)->schedIDX());
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

void Diff::reset(void)
{
}

////////////////////////////////////////////////////////////////////////////////

double Diff::rate(void) const
{
	// Pre-fetch some general info.
	CompDef * cdef = pTet->compdef();
    // Fetch the ligand as global index.
    uint gidx = pDiffDef->lig();
    // As local index.
    uint lidx = cdef->specG2L(gidx);
    
    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pTet->poolCount(lidx));
    assert(std::isnan(rate) == false);
    
    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

SchedIDXVec const & Diff::apply(State * s)
{
	// Pre-fetch some general info.
	CompDef * cdef = pTet->compdef();
    // Fetch the ligand as global index.
    uint gidx = def()->lig();
    // As local index.
    uint lidx = cdef->specG2L(gidx);
    
    // Apply local change.
    pTet->incPoolCount(lidx, -1);
    
    // Apply change in next voxel: select a direction.
    double sel = s->rng()->getUnfEE();
    if (sel < pCDFSelector[0])
    {
        // Direction 1.
    	Tet * next = pTet->nextTet(0);
        next->incPoolCount(lidx, 1);
        return pUpdVec[0];
    }
    else if (sel < pCDFSelector[1])
    {
        // Direction 2.
    	Tet * next = pTet->nextTet(1);
        next->incPoolCount(lidx, 1);
        return pUpdVec[1];
    }
    else if (sel < pCDFSelector[2])
    {
        // Direction 3.
    	Tet * next = pTet->nextTet(2);
        next->incPoolCount(lidx, 1);
        return pUpdVec[2];
    }
    else 
    {
        // Direction 4.
    	Tet * next = pTet->nextTet(3);
        next->incPoolCount(lidx, 1);
        return pUpdVec[3];
    }
    
    // This should never happen!
    assert(0);
    return pUpdVec[0];
}

////////////////////////////////////////////////////////////////////////////////

// END
