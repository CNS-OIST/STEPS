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
// $Id:func_core.cpp 64 2007-08-20 06:25:41Z stefan $
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
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>
#include <steps/sim/swiginf/func_core.hpp>
#include <steps/tetexact/solver_core/state.hpp>

NAMESPACE_ALIAS(steps::sim, ssim);
NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverName(void)
{
    return "tetexact";
}

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverDesc(void)
{
    return "Exact SSA in tetrahedral mesh";
}

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverAuthors(void)
{
    return "Stefan Wils";
}

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverEmail(void)
{
    return "wils@oist.jp";
}

////////////////////////////////////////////////////////////////////////////////

State * siNewState(void)
{
    return new State();
}

////////////////////////////////////////////////////////////////////////////////

void siDelState(State * s)
{
    delete s;
}

////////////////////////////////////////////////////////////////////////////////

void siBeginStateDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndStateDef(State * s)
{
    s->def()->setupFinal();
    s->setupState();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginVarDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndVarDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siNewSpec(State * s, char * name)
{
    ssim::SpecDef * spec = s->def()->createSpecDef(name);
    assert(spec != 0);
    return spec->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginReacDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndReacDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siNewReac(State * s, char * name, double kcst)
{
    ssim::ReacDef * reac = s->def()->createReacDef(name);
    assert(reac != 0);
    reac->setKcst(kcst);
    return reac->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siAddReacLHS(State * s, uint ridx, uint sidx)
{
    ssim::ReacDef * reac = s->def()->reac(ridx);
    assert(reac != 0);
    reac->incLHS(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddReacRHS(State * s, uint ridx, uint sidx)
{
    ssim::ReacDef * reac = s->def()->reac(ridx);
    assert(reac != 0);
    reac->incRHS(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginDiffDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndDiffDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siNewDiff(State * s, char * name, uint sidx, double dcst)
{
    ssim::DiffDef * diff = s->def()->createDiffDef(name);
    assert(diff != 0);
    diff->setDcst(dcst);
    diff->setLig(sidx);
    return diff->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginSReacDef(State * s)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siEndSReacDef(State * s)
{
    
}

////////////////////////////////////////////////////////////////////////////////

uint siNewSReac(State * s, char * name, double kcst, bool inside)
{
    ssim::SReacDef * sreac = s->def()->createSReacDef(name, inside);
    assert(sreac != 0);
    sreac->setKcst(kcst);
    return sreac->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacLHS_I(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incLHS_I(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacLHS_S(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incLHS_S(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacLHS_O(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incLHS_O(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacRHS_I(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incRHS_I(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacRHS_S(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incRHS_S(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddSReacRHS_O(State * s, uint ridx, uint sidx)
{
    ssim::SReacDef * sreac = s->def()->sreac(ridx);
    assert(sreac != 0);
    sreac->incRHS_O(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginCompDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

void siEndCompDef(State * s)
{
}

////////////////////////////////////////////////////////////////////////////////

uint siNewComp(State * s, char * name, double vol)
{
    ssim::CompDef * compdef = s->def()->createCompDef(name);
    assert(compdef != 0);
    compdef->setVol(vol);
    uint compdef_gidx = compdef->gidx();
    
    uint comp_idx = s->addComp(compdef);
    assert(compdef_gidx == comp_idx);
    
    return compdef_gidx;
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompSpec(State * s, uint cidx, uint sidx)
{
    ssim::CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addSpec(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompReac(State * s, uint cidx, uint ridx)
{
    ssim::CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addReac(ridx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompDiff(State * s, uint cidx, uint didx)
{
    ssim::CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addDiff(didx);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginPatchDef(State * s)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siEndPatchDef(State * s)
{
    
}

////////////////////////////////////////////////////////////////////////////////

uint siNewPatch
(
    State * s, 
    char * name, 
    double area, 
    uint cidx_in, 
    uint cidx_out
)
{
    ssim::CompDef * cdef_i = 0;
    if (cidx_in != 0xFFFF) cdef_i = s->def()->comp(cidx_in);
    ssim::CompDef * cdef_o = 0;
    if (cidx_out != 0xFFFF) cdef_o = s->def()->comp(cidx_out);
    
    ssim::PatchDef * pdef = s->def()->createPatchDef(name, cdef_i, cdef_o);
    assert(pdef != 0);
    pdef->setArea(area);
    uint patchdef_gidx = pdef->gidx();
    
    uint patch_idx = s->addPatch(pdef);
    assert(patchdef_gidx == patch_idx);
    
    return patchdef_gidx;
}

////////////////////////////////////////////////////////////////////////////////

void siAddPatchSpec(State * s, uint pidx, uint sidx)
{
    ssim::PatchDef * pdef = s->def()->patch(pidx);
    assert(pdef != 0);
    pdef->addSpec(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddPatchSReac(State * s, uint pidx, uint ridx)
{
    ssim::PatchDef * pdef = s->def()->patch(pidx);
    assert(pdef != 0);
    pdef->addSReac(ridx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetRNG(State * s, steps::rng::RNG * rng)
{
    s->setRNG(rng);
}

////////////////////////////////////////////////////////////////////////////////

void siReset(State * s)
{
    s->reset();
}

////////////////////////////////////////////////////////////////////////////////

void siRun(State * s, double endtime)
{
    s->run(endtime);
}

////////////////////////////////////////////////////////////////////////////////

double siGetTime(State * s)
{
    return s->time();
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompVol(State * s, uint cidx)
{
    assert(s != 0);
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    return comp->vol();
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompVol(State * s, uint cidx, double vol)
{
    // Not implemented!
}

////////////////////////////////////////////////////////////////////////////////

uint siGetCompCount(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    uint slidx = comp->def()->specG2L(sidx);
    if (slidx == ssim::LIDX_UNDEFINED) return 0;
    
    uint count = 0;
    TetPVecCI t_end = comp->endTet();
    for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
    {
        count += (*t)->pools()[slidx];
    }
    
    return count;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompCount(State * s, uint cidx, uint sidx, uint n)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    ssim::lidxT slidx = comp->def()->specG2L(sidx);
    if (slidx == ssim::LIDX_UNDEFINED) return;
    
    double totalvol = comp->vol();
    
    if (n >= comp->countTets())
    {
        TetPVecCI t_end = comp->endTet();
        for (TetPVecCI t = comp->bgnTet(); t != t_end; ++t)
        {
            Tet * tet = *t;
            double fract = static_cast<double>(n) * (tet->vol() / totalvol);
            uint n3 = static_cast<uint>(std::floor(fract));
            tet->pools()[slidx] = n3;
            n -= n3;
        }
    }
    
    while (n != 0)
    {
        Tet * tet = comp->pickTetByVol(s->rng()->getUnfIE());
        assert(tet != 0);
        tet->pools()[slidx] += 1;
        n--;
    }
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompMass(State * s, uint cidx, uint sidx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompMass(State * s, uint cidx, uint sidx, double m)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompConc(State * s, uint cidx, uint sidx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompConc(State * s, uint cidx, uint sidx, double c)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompClamped(State * s, uint cidx, uint sidx)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompClamped(State * s, uint cidx, uint sidx, bool buf)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacK(State * s, uint cidx, uint ridx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacK(State * s, uint cidx, uint ridx, double kf)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompReacActive(State * s, uint cidx, uint ridx)
{
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacActive(State * s, uint cidx, uint ridx, bool act)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompDiffD(State * s, uint cidx, uint didx)
{
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompDiffD(State * s, uint cidx, uint didx)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompDiffActive(State * s, uint cidx, uint didx)
{
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompDiffActive(State * s, uint cidx, uint didx, bool act)
{
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchArea(State * s, uint pidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchArea(State * s, uint pidx, double area)
{
    
}

////////////////////////////////////////////////////////////////////////////////

uint siGetPatchCount(State * s, uint pidx, uint sidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchCount(State * s, uint pidx, uint sidx, uint n)
{
    
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchMass(State * s, uint pidx, uint sidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchMass(State * s, uint pidx, uint sidx, double m)
{
}

////////////////////////////////////////////////////////////////////////////////

bool siGetPatchClamped(State * s, uint pidx, uint sidx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchClamped(State * s, uint pidx, uint sidx, bool buf)
{
    
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchSReacK(State * s, uint pidx, uint ridx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchSReacK(State * s, uint pidx, uint ridx, double kf)
{
    
}

////////////////////////////////////////////////////////////////////////////////

bool siGetPatchSReacActive(State * s, uint pidx, uint ridx)
{
    
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchSReacActive(State * s, uint pidx, uint ridx, bool a)
{
}

////////////////////////////////////////////////////////////////////////////////

// END
