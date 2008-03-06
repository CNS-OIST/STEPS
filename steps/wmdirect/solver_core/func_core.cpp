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
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/swiginf/func_core.hpp>
#include <steps/wmdirect/solver_core/comp.hpp>
#include <steps/wmdirect/solver_core/patch.hpp>
#include <steps/wmdirect/solver_core/reac.hpp>
#include <steps/wmdirect/solver_core/sched.hpp>
#include <steps/wmdirect/solver_core/sreac.hpp>
#include <steps/wmdirect/solver_core/state.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);
NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverName(void)
{
    return "wmdirect";
}

////////////////////////////////////////////////////////////////////////////////

char * siGetSolverDesc(void)
{
    return "SSA Direct Method in well-mixed conditions";
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

    // Create the actual compartments.
    ssim::CompDefPVecCI c_end = s->def()->endComp();
    for (ssim::CompDefPVecCI c = s->def()->bgnComp(); c != c_end; ++c)
    {
        uint compdef_gidx = (*c)->gidx();
        uint comp_idx = s->addComp(*c);
        assert(compdef_gidx == comp_idx);
    }
    
    // Create the actual patches.
    ssim::PatchDefPVecCI p_end = s->def()->endPatch();
    for (ssim::PatchDefPVecCI p = s->def()->bgnPatch(); p != p_end; ++p)
    {
        uint patchdef_gidx = (*p)->gidx();
        uint patch_idx = s->addPatch(*p);
        assert(patchdef_gidx == patch_idx);
    }
    
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
    return compdef->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompSpec(State * s, uint cidx, uint sidx)
{
    ssim::CompDef * compdef = s->def()->comp(cidx);
    assert(compdef != 0);
    compdef->addSpec(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompReac(State * s, uint cidx, uint ridx)
{
    ssim::CompDef * compdef = s->def()->comp(cidx);
    assert(compdef != 0);
    compdef->addReac(ridx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompDiff(State * s, uint cidx, uint didx)
{
    ssim::CompDef * compdef = s->def()->comp(cidx);
    assert(compdef != 0);
    compdef->addDiff(didx);
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
    
    return pdef->gidx();
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
    
    return comp->pools()[slidx];
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
    
    comp->pools()[slidx] = n;
    // It's cheaper to just recompute everything.
    s->sched()->reset();
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompAmount(State * s, uint cidx, uint sidx)
{
    double count = static_cast<double>(siGetCompCount(s, cidx, sidx));
    return count / smath::AVOGADRO; 
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompAmount(State * s, uint cidx, uint sidx, double m)
{
    assert(s != 0);
    assert(m >= 0.0);
    
    double m2 = m * smath::AVOGADRO;
    double m_int = std::floor(m2);
    double m_frc = m2 - m_int;
    uint c = static_cast<uint>(m_int);
    if (m_frc > 0.0)
    {
        double rand01 = s->rng()->getUnfIE();
        if (rand01 < m_frc) c++;
    }
    
    siSetCompCount(s, cidx, sidx, c);
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompConc(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    
    double count = static_cast<double>(siGetCompCount(s, cidx, sidx));
    double vol = comp->vol();
    return count / (1.0e3 * vol * smath::AVOGADRO); 
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompConc(State * s, uint cidx, uint sidx, double c)
{
    assert(s != 0);
    assert(c >= 0.0);

    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);    
    
    double c2 = c * (1.0e3 * comp->vol() * smath::AVOGADRO);
    double c_int = std::floor(c2);
    double c_frc = c2 - c_int;
    uint count = static_cast<uint>(c_int);
    if (c_frc > 0.0)
    {
        double rand01 = s->rng()->getUnfIE();
        if (rand01 < c_frc) count++;
    }
    
    siSetCompCount(s, cidx, sidx, count);
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompClamped(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    
    uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return false;
    
    return comp->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompClamped(State * s, uint cidx, uint sidx, bool buf)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    
    uint lsidx = comp->def()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return;
    
    comp->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacK(State * s, uint cidx, uint ridx)
{
    // Currently not implemented.
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacK(State * s, uint cidx, uint ridx, double kf)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompReacActive(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return false;
    
    return !(comp->reac(lridx)->inactive());
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacActive(State * s, uint cidx, uint ridx, bool act)
{
    assert(s != 0);
    assert(cidx < s->countComps());
    Comp * comp = s->comp(cidx);
    assert(comp != 0);
    
    uint lridx = comp->def()->reacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return;
    
    comp->reac(lridx)->setActive(act);
    
    // It's cheaper to just recompute everything.
    s->sched()->reset();
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
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    return patch->area();
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchArea(State * s, uint pidx, double area)
{
    // Not implemented!
}

////////////////////////////////////////////////////////////////////////////////

uint siGetPatchCount(State * s, uint pidx, uint sidx)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    uint slidx = patch->def()->specG2L(sidx);
    if (slidx == ssim::LIDX_UNDEFINED) return 0;
    
    return patch->pools()[slidx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchCount(State * s, uint pidx, uint sidx, uint n)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    uint slidx = patch->def()->specG2L(sidx);
    if (slidx == ssim::LIDX_UNDEFINED) return;
    
    patch->pools()[slidx] = n;
    // It's cheaper to just recompute everything.
    s->sched()->reset();
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchAmount(State * s, uint pidx, uint sidx)
{
    double count = static_cast<double>(siGetPatchCount(s, pidx, sidx));
    return count / smath::AVOGADRO; 
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchAmount(State * s, uint pidx, uint sidx, double m)
{
    assert(s != 0);
    assert(m >= 0.0);
    
    double m2 = m * smath::AVOGADRO;
    double m_int = std::floor(m2);
    double m_frc = m2 - m_int;
    uint c = static_cast<uint>(m_int);
    if (m_frc > 0.0)
    {
        double rand01 = s->rng()->getUnfIE();
        if (rand01 < m_frc) c++;
    }
    
    siSetPatchCount(s, pidx, sidx, c);
}

////////////////////////////////////////////////////////////////////////////////

bool siGetPatchClamped(State * s, uint pidx, uint sidx)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    
    uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return false;
    
    return patch->clamped(lsidx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchClamped(State * s, uint pidx, uint sidx, bool buf)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    
    uint lsidx = patch->def()->specG2L(sidx);
    if (lsidx == ssim::LIDX_UNDEFINED) return;
    
    patch->setClamped(lsidx, buf);
}

////////////////////////////////////////////////////////////////////////////////

double siGetPatchSReacK(State * s, uint pidx, uint ridx)
{
    // Currently not implemented.
    return 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchSReacK(State * s, uint pidx, uint ridx, double kf)
{
    // Currently not implemented.
}

////////////////////////////////////////////////////////////////////////////////

bool siGetPatchSReacActive(State * s, uint pidx, uint ridx)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    
    uint lridx = patch->def()->sreacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return false;
    
    return !(patch->sreac(lridx)->inactive());
}

////////////////////////////////////////////////////////////////////////////////

void siSetPatchSReacActive(State * s, uint pidx, uint ridx, bool a)
{
    assert(s != 0);
    assert(pidx < s->countPatches());
    Patch * patch = s->patch(pidx);
    assert(patch != 0);
    
    uint lridx = patch->def()->sreacG2L(ridx);
    if (lridx == ssim::LIDX_UNDEFINED) return;
    
    patch->sreac(lridx)->setActive(a);
    
    // It's cheaper to just recompute everything.
    s->sched()->reset();
}

////////////////////////////////////////////////////////////////////////////////

// END
