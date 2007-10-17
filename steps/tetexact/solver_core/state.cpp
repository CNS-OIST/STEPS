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
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

////////////////////////////////////////////////////////////////////////////////
    
State::State(void)
: pStateDef(0)
, pRNG(0)
, pSched(0)
, pTime(0.0)
, pTets()
{
    pStateDef = new StateDef();
    pSched = new Sched();
}

////////////////////////////////////////////////////////////////////////////////

State::~State(void)
{
    delete pStateDef;
}

////////////////////////////////////////////////////////////////////////////////

void State::setupState(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void State::setupTetmesh(void)
{
    // First we create all kinetic processes.
    for (std::vector<Tet*>::const_iterator i = pTets.begin(); 
        i != pTets.end(); ++i)
    {
        (*i)->setupKProcs(pSched);
    }
    // Next, we resolve all dependencies.
    for (std::vector<Tet*>::const_iterator i = pTets.begin(); 
            i != pTets.end(); ++i)
    {
        for (std::vector<Diff*>::const_iterator j = (*i)->diffBegin();
            j != (*i)->diffEnd(); ++j)
        {
            (*j)->setupDeps();
        }
    }
    sched()->build();
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
    std::for_each(pTets.begin(), pTets.end(), std::mem_fun(&Tet::reset));
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

uint State::addTet
(
    CompDef * cdef, double vol, 
    double a1, double a2, double a3, double a4,
    double d1, double d2, double d3, double d4
)
{
    Tet * t = new Tet(cdef, vol, a1, a2, a3, a4, d1, d2, d3, d4);
    uint tidx = pTets.size();
    pTets.push_back(t);
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
