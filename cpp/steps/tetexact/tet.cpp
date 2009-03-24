////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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
//
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/compdef.hpp>
#include <steps/solver/diffdef.hpp>
#include <steps/solver/reacdef.hpp>
#include <steps/tetexact/diff.hpp>
#include <steps/tetexact/reac.hpp>
#include <steps/tetexact/tet.hpp>
#include <steps/tetexact/tri.hpp>
#include <steps/tetexact/kproc.hpp>
#include <steps/tetexact/tetexact.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);
NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

stex::Tet::Tet
(
    solver::Compdef * cdef, double vol,
    double a0, double a1, double a2, double a3,
    double d0, double d1, double d2, double d3,
    int tet0, int tet1, int tet2, int tet3
)
: pCompdef(cdef)
, pVol(vol)
, pTets()
//, pTris()
, pNextTet()
, pNextTri()
, pAreas()
, pDist()
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
{
    assert(pCompdef != 0);
	assert (pVol > 0.0);
	assert (a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i=0; i <= 3; ++i)
    {
    	pNextTet[i] = 0;
    	pNextTri[i] = 0;
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

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    pKProcs.resize(compdef()->countDiffs() + compdef()->countReacs());

}

////////////////////////////////////////////////////////////////////////////////

stex::Tet::~Tet(void)
{
    // Delete species pool information.
    delete[] pPoolCount;
    delete[] pPoolFlags;

    // Delete diffusion rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTet(uint i, stex::Tet * t)
{
    if (t->compdef() != compdef())
    {
        pNextTet[i] = 0;
    }
    else
    {
        pNextTet[i] = t;
        if (pNextTri[i] != 0) std::cout << "WARNING: writing over nextTri index " << i;
        pNextTri[i] = 0;
    }
}


////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setNextTri(uint i, stex::Tri * t)
{
	if (pNextTet[i] != 0) std::cout << "WARNING: writing over nextTet index " << i;
    pNextTet[i] = 0;
    pNextTri[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setupKProcs(stex::Tetexact * tex)
{
    uint j = 0;
    // Create diffusion kproc's.
    uint ndiffs = compdef()->countDiffs();
    for (uint i = 0; i < ndiffs; ++i)
    {
        ssolver::Diffdef * ddef = compdef()->diffdef(i);
        stex::Diff * d = new stex::Diff(ddef, this);
        pKProcs[j++] = d;
        tex->addKProc(d);
    }

    // Create reaction kproc's.
    uint nreacs = compdef()->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        ssolver::Reacdef * rdef = compdef()->reacdef(i);
        stex::Reac * r = new stex::Reac(rdef, this);
        pKProcs[j++] = r;
        tex->addKProc(r);
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::reset(void)
{
    uint nspecs = compdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);
    std::for_each(pKProcs.begin(), pKProcs.end(),
    		std::mem_fun(&stex::KProc::reset));
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setCount(uint lidx, uint count)
{
	assert (lidx < compdef()->countSpecs());
	double oldcount = pPoolCount[lidx];
	double c = static_cast<double>(count);
	pPoolCount[lidx] = c;
	// Now update the count in this tet's comp
	double diff = c - oldcount;
	double newcount = (compdef()->pools()[lidx]) + diff;
	// Compdef method will do the checking on the double argument
	// (should be positive or zero!)
	compdef()->setCount(lidx, newcount);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tet::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

stex::Diff * stex::Tet::diff(uint lidx) const
{
    assert(lidx < compdef()->countDiffs());
    return dynamic_cast<stex::Diff*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////

stex::Reac * stex::Tet::reac(uint lidx) const
{
    assert(lidx < compdef()->countReacs());
    return dynamic_cast<stex::Reac*>(pKProcs[compdef()->countDiffs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

// END
