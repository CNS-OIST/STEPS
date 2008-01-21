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
// $Id:state.cpp 64 2007-08-20 06:25:41Z stefan $
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
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/tetexact/solver_core/comp.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/patch.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////
    
State::State(void)
: pStateDef(0)
, pRNG(0)
, pSched(0)
, pTime(0.0)
, pComps()
, pTets()
, pPatches()
, pTris()
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
    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI tet = pTets.begin(); tet != tet_e; ++tet) delete *tet;
    PatchPVecCI patch_e = pPatches.end();
    for (PatchPVecCI p = pPatches.begin(); p != patch_e; ++p) delete *p;
    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI tri = pTris.begin(); tri != tri_e; ++tri) delete *tri;
}

////////////////////////////////////////////////////////////////////////////////

void State::setupState(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void State::setupTetmesh(void)
{
    // First we create all kinetic processes.
    TetPVecCI tet_e = pTets.end();
    for (TetPVecCI tet = pTets.begin(); tet != tet_e; ++tet) 
    {
        (*tet)->setupKProcs(pSched);
    }
    TriPVecCI tri_e = pTris.end();
    for (TriPVecCI tri = pTris.begin(); tri != tri_e; ++tri)
    {
        (*tri)->setupKProcs(pSched);
    }
    
    // Next, we resolve all dependencies.
    for (TetPVecCI tet = pTets.begin(); tet != tet_e; ++tet)
    {
        KProcPVecCI kprocend = (*tet)->kprocEnd();
        for (KProcPVecCI k = (*tet)->kprocBegin(); k != kprocend; ++k)
        {
            (*k)->setupDeps();
        }
    }
    for (TriPVecCI tri = pTris.begin(); tri != tri_e; ++tri)
    {
        KProcPVecCI kprocend = (*tri)->kprocEnd();
        for (KProcPVecCI k = (*tri)->kprocBegin(); k != kprocend; ++k)
        {
            (*k)->setupDeps();
        }
    }
    
    sched()->build();
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
    std::for_each(pTets.begin(), pTets.end(), std::mem_fun(&Tet::reset));
    std::for_each(pTris.begin(), pTris.end(), std::mem_fun(&Tri::reset));
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

uint State::addComp(ssim::CompDef * cdef)
{
    Comp * comp = new Comp(cdef);
    assert(comp != 0);
    uint compidx = pComps.size();
    pComps.push_back(comp);
    return compidx;
}

////////////////////////////////////////////////////////////////////////////////

Comp * State::comp(uint idx) const
{
    assert(idx < pComps.size());
    return pComps[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint State::addTet
(
    Comp * comp, double vol, 
    double a1, double a2, double a3, double a4,
    double d1, double d2, double d3, double d4
)
{
    ssim::CompDef * compdef = comp->def(); 
    TetP t = new Tet(compdef, vol, a1, a2, a3, a4, d1, d2, d3, d4);
    uint tidx = pTets.size();
    pTets.push_back(t);
    comp->addTet(t);
    return tidx;
}

////////////////////////////////////////////////////////////////////////////////

uint State::addPatch(ssim::PatchDef * pdef)
{
    Patch * patch = new Patch(pdef);
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

uint State::addTri(Patch * patch, double area)
{
    ssim::PatchDef * patchdef = patch->def();
    TriP t = new Tri(patchdef, area);
    uint tidx = pTris.size();
    pTris.push_back(t);
    patch->addTri(t);
    return tidx;
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
