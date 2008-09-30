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
#include <steps/sim/shared/patchdef.hpp>														
#include <steps/wmrk4/solver_core/patch.hpp>
#include <steps/wmrk4/solver_core/comp.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

Patch::Patch(ssim::PatchDef * patchdef, Comp * icomp, Comp * ocomp)
: pPatchDef(patchdef)
, pPoolCount(0)																			
, pPoolFlags(0)
, pSReacFlags(0)
, pIComp(icomp)
, pOComp(ocomp)
{
    assert(pPatchDef != 0);
    uint nspecs = def()->countSpecs();
	uint nreacs = def()->countSReacs();
    pPoolCount = new double[nspecs];														
    pPoolFlags = new uint[nspecs];
	pSReacFlags = new uint[nreacs];
    std::fill_n(pPoolCount, nspecs, 0.0);
    std::fill_n(pPoolFlags, nspecs, 0);
	std::fill_n(pSReacFlags, nreacs, 0);
    if (iComp() != 0) iComp()->addIPatch(this);
    if (oComp() != 0) oComp()->addOPatch(this);
}

////////////////////////////////////////////////////////////////////////////////

Patch::~Patch(void)
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
	delete[] pSReacFlags;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::reset(void)
{
    uint nspecs = def()->countSpecs();
	uint nreacs = def()->countSReacs();
    std::fill_n(pPoolCount, nspecs, 0.0);															
    std::fill_n(pPoolFlags, nspecs, 0);
	std::fill_n(pSReacFlags, nreacs, 0);
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

void Patch::setActive(uint lidx, bool active)
{
	if (active == true) pSReacFlags[lidx] &= ~INACTIVATED;
	else pSReacFlags[lidx] |= INACTIVATED;
}
////////////////////////////////////////////////////////////////////////////////


// END


