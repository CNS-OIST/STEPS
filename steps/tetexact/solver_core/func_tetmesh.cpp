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
// $Id:func_tetmesh.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/swiginf/func_ssa.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

////////////////////////////////////////////////////////////////////////////////

void siBeginTetmeshDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTetmeshDef(State *s)
{
    s->setupTetmesh();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginTetDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTetDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siNewTet(State * s, uint cidx, double vol, 
    double a1, double a2, double a3, double a4,
    double d1, double d2, double d3, double d4)
{
    CompDef * cdef = s->def()->comp(cidx);
    return s->addTet(cdef, vol, a1, a2, a3, a4, d1, d2, d3, d4);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginConnectDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndConnectDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTetTet(State * s, uint side, uint tidx1, uint tidx2)
{
    Tet * t1 = s->tet(tidx1);
    Tet * t2 = s->tet(tidx2);
    t1->setNextTet(side, t2);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTetTriInside(State * s, uint side, uint tetidx, uint triidx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTetTriOutside(State * s, uint side, uint tetidx, uint triidx)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetVol(State * s, uint tidx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetVol(State * s, uint tidx, double vol)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siGetTetCount(State * s, uint tidx, uint sidx)
{
    Tet * tet = s->tet(tidx);
    // TODO: error stuff
    if (tet == 0) return 0;
    uint l_sidx = tet->compdef()->specG2L(sidx);
    if (l_sidx == 0xFFFF) return 0;
    return tet->poolCount(l_sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetCount(State * s, uint tidx, uint sidx, uint n)
{
    Tet * tet = s->tet(tidx);
    // TODO: error stuff
    if (tet == 0) return;
    
    // Apply the change.
    uint l_sidx = tet->compdef()->specG2L(sidx);
    if (l_sidx == 0xFFFF) return;
    tet->setPoolCount(l_sidx, n);
    
    // Make updates to the schedule.
    // DEBUG: 04-Sep-2007
    CompUpd * cupd = tet->compdef()->updateSpec(l_sidx);
    SchedIDXVec updvec;
    
    // Loop over diffusions.
    std::vector<uint>::const_iterator diff_end = cupd->endLDiffs();
    for (std::vector<uint>::const_iterator diff = cupd->beginLDiffs();
        diff != diff_end; ++diff)
    {
        updvec.push_back(tet->diff(*diff)->schedIDX());
    }
    
    // Loop over reactions (not yet).
    //
    
    // Send the list of kprocs that need to be updated to the schedule.
    s->sched()->update(updvec);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetMass(State * s, uint tidx, uint sidx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetMass(State * s, uint tidx, uint sidx, double m)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetConc(State * s, uint tidx, uint sidx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetConc(State * s, uint tidx, uint sidx, double c)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetClamped(State * s, uint tidx, uint sidx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetClamped(State * s, uint tidx, uint sidx, bool buf)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetReacK(State * s, uint tidx, uint ridx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetReacK(State * s, uint tidx, uint ridx, double kf)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetReacActive(State * s, uint tidx, uint ridx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetReacActive(State * s, uint tidx, uint ridx, bool act)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetDiffD(State * s, uint tidx, uint didx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siSetTetDiffD(State * s, uint tidx, uint didx)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTetDiffActive(State * s, uint tidx, uint didx)
{
}

////////////////////////////////////////////////////////////////////////////////

void siGetTetDiffActive(State * s, uint tidx, uint didx, bool act)
{
}

////////////////////////////////////////////////////////////////////////////////

// END