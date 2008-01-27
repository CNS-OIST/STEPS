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
#include <steps/wmdirect/solver_core/comp.hpp>
#include <steps/wmdirect/solver_core/kproc.hpp>
#include <steps/wmdirect/solver_core/patch.hpp>
#include <steps/wmdirect/solver_core/sched.hpp>
#include <steps/wmdirect/solver_core/state.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

State::State(void)
: pStateDef(0)
, pRNG(0)
, pTime(0.0)
, pNSteps(0) 
, pSched(0)
, pComps()
, pCompMap()   
, pPatches()
{
    pStateDef = new ssim::StateDef();
    pSched = new Sched();
}

////////////////////////////////////////////////////////////////////////////////

State::~State(void)
{
    delete pStateDef;
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI c = pComps.begin(); c != comp_e; ++c) delete *c;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;
}

////////////////////////////////////////////////////////////////////////////////

void State::setupState(void)
{
    // First we create all kinetic processes.
    CompPVecCI comp_e = pComps.end();
    for (CompPVecCI comp = pComps.begin(); comp != comp_e; ++comp) 
    {
        (*comp)->setupKProcs(pSched);
    }
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI patch = pPatches.begin(); patch != patch_e; ++patch)
    {
        (*patch)->setupKProcs(pSched);
    }
    
    // Next, we resolve all dependencies.
    for (CompPVecCI comp = pComps.begin(); comp != comp_e; ++comp)
    {
        KProcPVecCI kprocend = (*comp)->kprocEnd();
        for (KProcPVecCI k = (*comp)->kprocBegin(); k != kprocend; ++k)
        {
            (*k)->setupDeps();
        }
    }
    for (PatchPVecCI patch = pPatches.begin(); patch != patch_e; ++patch)
    {
        KProcPVecCI kprocend = (*patch)->kprocEnd();
        for (KProcPVecCI k = (*patch)->kprocBegin(); k != kprocend; ++k)
        {
            (*k)->setupDeps();
        }
    }
    
    sched()->build();
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
    std::for_each(pComps.begin(), pComps.end(), std::mem_fun(&Comp::reset));
    std::for_each(pPatches.begin(), pPatches.end(), std::mem_fun(&Patch::reset));
    pSched->reset();
    resetTime();
    resetNSteps();
}

////////////////////////////////////////////////////////////////////////////////

void State::step(void)
{
    KProc * kp = sched()->getNext(this);
    if (kp == 0) return;
    double a0 = sched()->getA0();
    if (a0 == 0.0) return;
    double dt = rng()->getExp(a0);
    executeStep(kp, dt);
}

////////////////////////////////////////////////////////////////////////////////

void State::run(double maxt)
{
    assert(maxt >= 0.0);
    while (time() < maxt)
    {
        KProc * kp = sched()->getNext(this);
        if (kp == 0) break;
        double a0 = sched()->getA0();
        if (a0 == 0.0) break;
        double dt = rng()->getExp(a0);
        if ((time() + dt) > maxt) break;
        executeStep(kp, dt);
    }
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

void State::executeStep(KProc * kp, double dt)
{
    SchedIDXVec const & upd = kp->apply(this);
    sched()->update(upd);
    incTime(dt);
    incNSteps(1);
}

////////////////////////////////////////////////////////////////////////////////

// END
