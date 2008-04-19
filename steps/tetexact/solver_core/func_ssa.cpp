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
// $Id:func_ssa.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <cmath>
#include <iostream>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>
#include <steps/sim/swiginf/func_ssa.hpp>
#include <steps/tetexact/solver_core/comp.hpp>
#include <steps/tetexact/solver_core/diff.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/patch.hpp>
#include <steps/tetexact/solver_core/reac.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/sreac.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

double siStep(State * s)
{
    s->step();
    return s->time();
}

////////////////////////////////////////////////////////////////////////////////

uint siGetNSteps(State * s)
{
    return s->nsteps();
}

////////////////////////////////////////////////////////////////////////////////

double siGetA0(State * s)
{
    return s->sched()->getA0();
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacC(State * s, uint cidx, uint ridx)
{
    Comp * c = s->comp(cidx);
    CompDef * cdef = c->def();
    lidxT lridx = cdef->reacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0.0;
    
    TetPVecCI t_bgn = c->bgnTet();
    TetPVecCI t_end = c->endTet();
    if (t_bgn == t_end) return 0.0;
    double c2 = 0.0;
    double v = 0.0;
    for (TetPVecCI t = t_bgn; t != t_end; ++t)
    {
        double v2 = (*t)->vol();
        Reac * reac = (*t)->reac(lridx);
        c2 += reac->c() * v2;
        v += v2;
    }
    
    assert(v > 0.0);
    return c2 / v; 
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacH(State * s, uint cidx, uint ridx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacA(State * s, uint cidx, uint ridx)
{
    Comp * c = s->comp(cidx);
    CompDef * cdef = c->def();
    lidxT lridx = cdef->reacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0.0;
    
    TetPVecCI t_bgn = c->bgnTet();
    TetPVecCI t_end = c->endTet();
    if (t_bgn == t_end) return 0.0;
    double a = 0.0;
    for (TetPVecCI t = t_bgn; t != t_end; ++t)
    {
        Reac * reac = (*t)->reac(lridx);
        a += reac->rate();
    }
    
    return a; 
}

////////////////////////////////////////////////////////////////////////////////

uint siGetCompReacExtent(State * s, uint cidx, uint ridx)
{
    Comp * c = s->comp(cidx);
    CompDef * cdef = c->def();
    lidxT lridx = cdef->reacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0;
        
    TetPVecCI t_bgn = c->bgnTet();
    TetPVecCI t_end = c->endTet();
    if (t_bgn == t_end) return 0;
    uint x = 0;
    for (TetPVecCI t = t_bgn; t != t_end; ++t)
    {
        Reac * reac = (*t)->reac(lridx);
        x += reac->getExtent();
    }
        
    return x; 
}

////////////////////////////////////////////////////////////////////////////////

void siResetCompReacExtent(State * s, uint cidx, uint ridx)
{
    Comp * c = s->comp(cidx);
    CompDef * cdef = c->def();
    lidxT lridx = cdef->reacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return;
        
    TetPVecCI t_bgn = c->bgnTet();
    TetPVecCI t_end = c->endTet();
    if (t_bgn == t_end) return;
    for (TetPVecCI t = t_bgn; t != t_end; ++t)
    {
        Reac * reac = (*t)->reac(lridx);
        reac->resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchSReacC(State * s, uint pidx, uint ridx)
{
    Patch * p = s->patch(pidx);
    PatchDef * pdef = p->def();
    lidxT lridx = pdef->sreacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0.0;
        
    TriPVecCI t_bgn = p->bgnTri();
    TriPVecCI t_end = p->endTri();
    if (t_bgn == t_end) return 0.0;
    double c = 0.0;
    double a = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        double a2 = (*t)->area();
        SReac * sreac = (*t)->sreac(lridx);
        c += sreac->c() * a2;
        a += a2;
    }
    
    assert(a > 0.0);
    return c / a; 
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchSReacH(State * s, uint pidx, uint ridx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchSReacA(State * s, uint pidx, uint ridx)
{
    Patch * p = s->patch(pidx);
    PatchDef * pdef = p->def();
    lidxT lridx = pdef->sreacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0.0;
        
    TriPVecCI t_bgn = p->bgnTri();
    TriPVecCI t_end = p->endTri();
    if (t_bgn == t_end) return 0.0;
    double a = 0.0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        SReac * sreac = (*t)->sreac(lridx);
        a += sreac->rate();
    }
    
    return a; 
}

////////////////////////////////////////////////////////////////////////////////

uint siGetPatchSReacExtent(State * s, uint pidx, uint ridx)
{
    Patch * p = s->patch(pidx);
    PatchDef * pdef = p->def();
    lidxT lridx = pdef->sreacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return 0;
        
    TriPVecCI t_bgn = p->bgnTri();
    TriPVecCI t_end = p->endTri();
    if (t_bgn == t_end) return 0;
    uint x = 0;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        SReac * sreac = (*t)->sreac(lridx);
        x += sreac->getExtent();
    }
    
    return x;
}

////////////////////////////////////////////////////////////////////////////////

void siResetPatchSReacExtent(State * s, uint pidx, uint ridx)
{
    Patch * p = s->patch(pidx);
    PatchDef * pdef = p->def();
    lidxT lridx = pdef->sreacG2L(ridx);
    if (lridx == LIDX_UNDEFINED) return;
        
    TriPVecCI t_bgn = p->bgnTri();
    TriPVecCI t_end = p->endTri();
    if (t_bgn == t_end) return;
    for (TriPVecCI t = t_bgn; t != t_end; ++t)
    {
        SReac * sreac = (*t)->sreac(lridx);
        sreac->resetExtent();
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
