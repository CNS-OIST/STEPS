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
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/swiginf/func_tetmesh.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/reac.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);

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
    assert(s != 0);
    CompPVecCI c_end = s->endComp();
    for (CompPVecCI c = s->bgnComp(); c != c_end; ++c)
    {
        (*c)->computeVol();
    }
}

////////////////////////////////////////////////////////////////////////////////

uint siNewTet(State * s, uint cidx, double vol, 
    double a1, double a2, double a3, double a4,
    double d1, double d2, double d3, double d4)
{
    assert(s != 0);
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    return s->addTet(comp, vol, a1, a2, a3, a4, d1, d2, d3, d4);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginTriDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndTriDef(State * s)
{
    assert(s != 0);
    PatchPVecCI p_end = s->endPatch();
    for (PatchPVecCI p = s->bgnPatch(); p != p_end; ++p)
    {
        (*p)->computeArea();
    }
}

////////////////////////////////////////////////////////////////////////////////

uint siNewTri(State * s, uint pidx, double area)
{
    assert(s != 0);
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    assert(area >= 0.0);
    return s->addTri(patch, area);
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

void siConnectTetTri(State * s, uint side, uint tetidx, uint triidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tet->setNextTri(side, tri);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTriTetInner(State * s, uint triidx, uint tetidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tri->setInnerTet(tet);
}

////////////////////////////////////////////////////////////////////////////////

void siConnectTriTetOuter(State * s, uint triidx, uint tetidx)
{
    TetP tet = s->tet(tetidx);
    assert(tet != 0);
    TriP tri = s->tri(triidx);
    assert(tri != 0);
    tri->setOuterTet(tet);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTetVol(State * s, uint tidx)
{
    TetP tet = s->tet(tidx);
    assert(tet != 0);
    return tet->vol();
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
    return tet->pools()[l_sidx];
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
    tet->pools()[l_sidx] = n;
    
    // Send the list of kprocs that need to be updated to the schedule.
    s->sched()->updateSpec(tet, l_sidx);
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

double siGetTriArea(State * s, uint tidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriArea(State * s, uint tidx, double area)
{
    
}

////////////////////////////////////////////////////////////////////////////////

uint siGetTriCount(State * s, uint tidx, uint sidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriCount(State * s, uint tidx, uint sidx, uint n)
{
    
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTriClamped(State * s, uint tidx, uint sidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriClamped(State * s, uint tidx, uint sidx, bool buf)
{
    
}

////////////////////////////////////////////////////////////////////////////////

double siGetTriSReacK(State * s, uint tidx, uint ridx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriSReacK(State * s, uint tidx, uint ridx, double kf)
{
    
}

////////////////////////////////////////////////////////////////////////////////

bool siGetTriSReacActive(State * s, uint tidx, uint ridx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetTriSReacActive(State * s, uint tidx, uint ridx, bool act)
{
    
}

////////////////////////////////////////////////////////////////////////////////

// END
