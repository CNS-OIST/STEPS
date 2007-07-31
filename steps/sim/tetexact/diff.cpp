////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
// 
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <cmath>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/tetexact/diff.hpp>
#include <steps/sim/tetexact/kproc.hpp>
#include <steps/sim/tetexact/sched.hpp>
#include <steps/sim/tetexact/state.hpp>
#include <steps/sim/tetexact/tet.hpp>

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
	for (uint i = 0; i < 4; ++i)
	{
		Tet * next = pTet->nextTet(i);
		if (next == 0) continue;
		// Add yourself to the list of stuff.
		pUpdVec[i].push_back(schedIDX());
		pUpdVec[i].push_back(next->TMPDIFFIDX());
	}
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
    uint gidx = pDiffDef->lig();
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
