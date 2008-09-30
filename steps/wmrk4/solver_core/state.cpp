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
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/wmrk4/solver_core/comp.hpp>
#include <steps/wmrk4/solver_core/patch.hpp>
#include <steps/wmrk4/solver_core/state.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

State::State(void)
: pStateDef(0)
, pRNG(0)
, pTime(0.0)
, pWmrk4(0)
, pComps()
, pCompMap()   
, pPatches()
, pdt(0.0)															
{
    pStateDef = new ssim::StateDef();	
}

////////////////////////////////////////////////////////////////////////////////

State::~State(void)
{
    delete pStateDef;
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;
	delete pWmrk4;
}

////////////////////////////////////////////////////////////////////////////////

void State::setupState(void)																	
{	
	pWmrk4 = new Wmrk4(this);			
	assert(pWmrk4 != 0);
    wmrk4()->setup();
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
    std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
    std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));
	// recompute flags and values vectors in Wmrk4 object
	wmrk4()->refill();																			
    resetTime();
}

////////////////////////////////////////////////////////////////////////////////

void State::step(void)  
{
	assert (pdt > 0.0);
    wmrk4()->rksteps(time(), time() + pdt);
}

////////////////////////////////////////////////////////////////////////////////

void State::run(double maxt)
{
    assert(maxt >= 0.0);
    wmrk4()->rksteps(time(), maxt);
    setTime(maxt);
}

////////////////////////////////////////////////////////////////////////////////

uint State::addComp(steps::sim::CompDef * cdef)
{
    Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    pCompMap[cdef] = comp;
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

Comp * State::comp(uint idx) const 
{
    assert(idx < pComps.size());
    return pComps[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint State::addPatch(steps::sim::PatchDef * pdef)
{
    Comp * icomp = 0;
    Comp * ocomp = 0;
    if (pdef->icompdef()) icomp = pCompMap[pdef->icompdef()];
    if (pdef->ocompdef()) ocomp = pCompMap[pdef->ocompdef()];
    Patch * patch = new Patch(pdef, icomp, ocomp);
    assert(patch != 0);
    uint patchidx = pPatches.size();
    pPatches.push_back(patch);
    return patchidx;
}

////////////////////////////////////////////////////////////////////////////////

Patch * State::patch(uint idx) const
{
    assert(idx < pPatches.size());
    return pPatches[idx];
}

////////////////////////////////////////////////////////////////////////////////

Wmrk4 * State::wmrk4(void) const
{
	assert (pWmrk4 != 0);
	return pWmrk4;
}
////////////////////////////////////////////////////////////////////////////////


// END


