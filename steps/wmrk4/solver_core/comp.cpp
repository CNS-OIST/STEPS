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
//  
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <algorithm>
#include <cassert>
#include <functional>
#include <iostream>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/wmrk4/solver_core/comp.hpp>														

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Comp::Comp(ssim::CompDef * compdef)
: pCompDef(compdef)
, pPoolCount(0)																					
, pPoolFlags(0)
, pReacFlags(0)
, pIPatches()
, pOPatches()
{
    assert(pCompDef != 0);
    uint nspecs = def()->countSpecs();
	uint nreacs = def()->countReacs();
    pPoolCount = new double[nspecs];																
    pPoolFlags = new uint[nspecs];
	pReacFlags = new uint[nreacs];
    std::fill_n(pPoolCount, nspecs, 0.0);															
    std::fill_n(pPoolFlags, nspecs, 0);
	std::fill_n(pReacFlags, nreacs, 0);
}

////////////////////////////////////////////////////////////////////////////////

Comp::~Comp(void)
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
	delete[] pReacFlags;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::reset(void)
{
    uint nspecs = def()->countSpecs();
	uint nreacs = def()->countReacs();
    std::fill_n(pPoolCount, nspecs, 0.0);															
    std::fill_n(pPoolFlags, nspecs, 0);
	std::fill_n(pReacFlags, nreacs, 0);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::setActive(uint lidx, bool active)
{
	if (active == true) pReacFlags[lidx] &= ~INACTIVATED;
	else pReacFlags[lidx] |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addIPatch(Patch * p)
{
    assert(std::find(pIPatches.begin(), pIPatches.end(), p) == pIPatches.end());
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void Comp::addOPatch(Patch * p)
{
    assert(std::find(pOPatches.begin(), pOPatches.end(), p) == pOPatches.end());
    pOPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////


// END

