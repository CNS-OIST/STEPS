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
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <cassert>
#include <algorithm>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/tetexact/diff.hpp>
#include <steps/sim/tetexact/sched.hpp>
#include <steps/sim/tetexact/tet.hpp>

////////////////////////////////////////////////////////////////////////////////

Tet::Tet
(
	CompDef * cdef, double vol, 
	double a0, double a1, double a2, double a3, 
	double d0, double d1, double d2, double d3
)
{
	// Copy all this crap.
	pCompDef = cdef;
	pNextTet[0] = 0;
	pNextTet[1] = 0;
	pNextTet[2] = 0;
	pNextTet[3] = 0;
	// Tetrahedral volumes.
	assert(vol >= 0.0);
	pVol = vol;
	// Areas of boundary triangles.
	assert(a0 >= 0.0);
	pAreas[0] = a0;
	assert(a1 >= 0.0);
	pAreas[1] = a1;
	assert(a2 >= 0.0);
	pAreas[2] = a2;
	assert(a3 >= 0.0);
	pAreas[3] = a3;
	// Distances to neighbouring tetrahedrons.
	assert(d0 >= 0.0);
	pDist[0] = d0;
	assert(d1 >= 0.0);
	pDist[1] = d1;
	assert(d2 >= 0.0);
	pDist[2] = d2;
	assert(d3 >= 0.0);
	pDist[3] = d3;
	
	// Based on compartment definition, build other structures.
	uint nspecs = compdef()->countSpecs();
	pPoolCount = new uint[nspecs];
	pPoolFlags = new uint[nspecs];
	std::fill_n(pPoolCount, nspecs, 0);
	std::fill_n(pPoolFlags, nspecs, 0);
	
	uint ndiffs = compdef()->countDiffs();
	for (uint i = 0; i < ndiffs; ++i)
		pDiffs.push_back(0);
}

////////////////////////////////////////////////////////////////////////////////

Tet::~Tet(void)
{
	// Delete species pool information.
	delete[] pPoolCount;
	delete[] pPoolFlags;
	
	// Delete diffusion rules.
	uint ndiffs = compdef()->countDiffs();
	for (uint i = 0; i < ndiffs; ++i)
		delete pDiffs[i];
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupKProcs(Sched * sched)
{
	// Create diffusion kproc's.
	uint ndiffs = compdef()->countDiffs();
	for (uint i = 0; i < ndiffs; ++i)
	{
		DiffDef * ddef = compdef()->diff(i);
		Diff * d = pDiffs[i] = new Diff(ddef, this);
		sched->addKProc(d);
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setupDeps(void)
{
	uint ndiffs = compdef()->countDiffs();
	for (uint i = 0; i < ndiffs; ++i)
	{
		pDiffs[i]->setupDeps();
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tet::reset(void)
{
	uint nspecs = compdef()->countSpecs();
	std::fill_n(pPoolCount, nspecs, 0);
	std::fill_n(pPoolFlags, nspecs, 0);
	
	uint ndiffs = compdef()->countDiffs();
	for (uint i = 0; i < ndiffs; ++i)
	{
		pDiffs[i]->reset();
	}
}

////////////////////////////////////////////////////////////////////////////////

void Tet::setNextTet(uint i, Tet * t)
{
	pNextTet[i] = t;
}
	
////////////////////////////////////////////////////////////////////////////////

// END
