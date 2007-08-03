////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
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
#include <steps/sim/wmdirect/state.hpp>

////////////////////////////////////////////////////////////////////////////////

// STEPS library.
NAMESPACE_ALIAS(steps::math, smath);

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
    assert(s != 0);
    // Perhaps in the future, this might also work from other modes,
    // to allow changes in the state to occur during simulation.
    assert(s->def()->mode() == StateDef::NAIVE_MODE);
    
    s->def()->setMode(StateDef::SETUP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siEndStateDef(State * s)
{
    assert(s != 0);
    StateDef * sdef = s->def();
    assert(s->def()->mode() == StateDef::SETUP_MODE);
    
    // Create the actual state, based on its definition.
    s->def()->setupFinal();
    s->setupState();
    
    // Finishing the setup mode automatically goes into the READY mode;
    // for now at least!
    s->def()->setMode(StateDef::READY_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginVarDef(State * s)
{
    assert(s != 0);
    // The SETUP_VAR mode is a submode of mode SETUP, and can only be accessed
    // from this mode.
    assert(s->def()->mode() == StateDef::SETUP_MODE);
    
    s->def()->setMode(StateDef::SETUP_VAR_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siEndVarDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_VAR_MODE);
    
    s->def()->setMode(StateDef::SETUP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

uint siNewSpec(State * s, char * name)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_VAR_MODE);
    assert(name != 0);
    
    SpecDef * spec = s->def()->createSpecDef(name);
    assert(spec != 0);
    return spec->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginReacDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_MODE);
    
    s->def()->setMode(StateDef::SETUP_REAC_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siEndReacDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_REAC_MODE);
    
    s->def()->setMode(StateDef::SETUP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

uint siNewReac(State * s, char * name, double kcst)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_REAC_MODE);
    assert(name != 0);
    
    ReacDef * reac = s->def()->createReacDef(name);
    assert(reac != 0);
    reac->setKcst(kcst);
    return reac->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siAddReacLHS(State * s, uint ridx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_REAC_MODE);
    assert(s->def()->isValidReac(ridx));
    assert(s->def()->isValidSpec(sidx));
    
    ReacDef * reac = s->def()->reac(ridx);
    assert(reac != 0);
    reac->incLHS(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddReacRHS(State * s, uint ridx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_REAC_MODE);
    assert(s->def()->isValidReac(ridx));
    assert(s->def()->isValidSpec(sidx));
    
    ReacDef * reac = s->def()->reac(ridx);
    assert(reac != 0);
    reac->incRHS(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siBeginDiffDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_MODE);
    
    s->def()->setMode(StateDef::SETUP_DIFF_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siEndDiffDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_DIFF_MODE);
    
    s->def()->setMode(StateDef::SETUP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

uint siNewDiff(State * s, char * name, uint sidx, double dcst)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_DIFF_MODE);
    assert(name != 0);
    
    DiffDef * diff = s->def()->createDiffDef(name);
    assert(diff != 0);
    diff->setDcst(dcst);
    diff->setLig(sidx);
    return diff->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siBeginCompDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_MODE);
    
    s->def()->setMode(StateDef::SETUP_COMP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

void siEndCompDef(State * s)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_COMP_MODE);
    
    s->def()->setMode(StateDef::SETUP_MODE);
}

////////////////////////////////////////////////////////////////////////////////

uint siNewComp(State * s, char * name, double vol)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_COMP_MODE);
    assert(name != 0);
    
    CompDef * comp = s->def()->createCompDef(name);
    assert(comp != 0);
    comp->setVol(vol);
    return comp->gidx();
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompSpec(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_COMP_MODE);
    assert(s->def()->isValidSpec(sidx));
    
    CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addSpec(sidx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompReac(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_COMP_MODE);
    assert(s->def()->isValidReac(ridx));
    
    CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addReac(ridx);
}

////////////////////////////////////////////////////////////////////////////////

void siAddCompDiff(State * s, uint cidx, uint didx)
{
    assert(s != 0);
    assert(s->def()->mode() == StateDef::SETUP_COMP_MODE);
    assert(s->def()->isValidDiff(didx));
    
    CompDef * comp = s->def()->comp(cidx);
    assert(comp != 0);
    comp->addDiff(didx);
}

////////////////////////////////////////////////////////////////////////////////

void siSetRNG(State * s, steps::rng::RNG * rng)
{
    assert(s != 0);
    assert(rng != 0);
    
    s->fRNG = rng;
}

////////////////////////////////////////////////////////////////////////////////

void siReset(State * s)
{
    assert(s != 0);
    s->reset();
}

////////////////////////////////////////////////////////////////////////////////

void siRun(State * s, double endtime)
{
    // Prefetch the number of compartments.
    uint ncomps = s->def()->countComps();
    assert(s != 0);
    while (s->fTime < endtime)
    {
        // Compute time of next reaction.
        double a0 = s->computeZeroProp();
        assert(a0 >= 0.0);
        if (a0 == 0.0)
        {
            s->fTime = endtime;
            return;
        }
        
        // Find out when the next reaction is going to occur.
        double tnext = s->fRNG->getExp(a0);

        // If it's too late, don't execute it.
        if ((s->fTime + tnext) > endtime)
        {
            s->fTime = endtime;
            return;
        }

        // Find out which reaction it is by finding the right interval.
        double ival = a0 * s->fRNG->getUnfIE();
        
        double accum = 0.0;
        uint ncomps = s->def()->countComps();
        uint nreacs = s->def()->comp(0)->countReacs();
        uint comp = 0;
        uint reac = 0;
        while (ival >= accum + s->fReacProps[comp][reac])
        {
            // Increase the accumulator.
            accum += s->fReacProps[comp][reac];
            // Go to the next reaction in the current compartment.
            ++reac;
            // If it was the last reaction, go to the next compartment.
            if (reac == nreacs)
            {
                ++comp;
                // Should never happen...?
                if (comp == ncomps) assert(0);
                reac = 0;
                nreacs = s->def()->comp(comp)->countReacs();
            }
        }
        
        // Execute the selected reaction.
        CompDef * cdef = s->def()->comp(comp);
        int * upd_vec = cdef->reacSpecUpds(reac);
        uint npools = cdef->countSpecs();
        for (uint pool = 0; pool < npools; ++pool)
        {
            if ((s->fPoolFlags[comp][pool] & State::CLAMPED_POOLFLAG) != 0) 
            {
                continue;
            }
            int nval = 
                static_cast<int>(s->fPoolCount[comp][pool]) + upd_vec[pool];
            assert(nval >= 0);
            s->fPoolCount[comp][pool] = static_cast<uint>(nval);
        }
        // Update extent.
        s->fReacExtents[comp][reac] += 1;
        
        // Update propensities for the selected reaction.
        CompUpd * cupd = cdef->updateReac(reac);
        for (std::vector<uint>::const_iterator i = cupd->beginLReacs();
            i != cupd->endLReacs(); ++i)
        {
            s->computeReacProp(comp, *i);
        }
        
        // Update time and statistics.
        s->fTime += tnext;
        s->fNSteps++;
    }
}

////////////////////////////////////////////////////////////////////////////////

double siGetTime(State * s)
{
    assert(s != 0);
    return s->fTime;
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompVol(State * s, uint cidx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    return s->fCompVols[cidx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompVol(State * s, uint cidx, double vol)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(vol >= 0.0);
    
    s->fCompVols[cidx] = vol;
    // No need to change local pools, because they're stored as numbers
    // of molecules. We do need to recompute all propensities though.
    uint nreacs = s->def()->comp(cidx)->countReacs();
    for (uint i = 0; i < nreacs; ++i)
    {
        // Recompute scaled reaction constants.
        s->computeCcst(cidx, i);
        // And of course also the corresponding propensity.
        s->computeReacProp(cidx, i);
    }
}

////////////////////////////////////////////////////////////////////////////////

uint siGetCompCount(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));
    
    uint l_sidx = s->def()->comp(cidx)->specG2L(sidx);
    if (l_sidx == 0xFFFF) return 0;
    return s->fPoolCount[cidx][l_sidx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompCount(State * s, uint cidx, uint sidx, uint n)
{

    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));

    CompDef * cdef = s->def()->comp(cidx);
    uint l_sidx = cdef->specG2L(sidx);
    if (l_sidx == 0xFFFF) return;
    s->fPoolCount[cidx][l_sidx] = n;

    // Recompute propensities.
    CompUpd * cupd = cdef->updateSpec(l_sidx);
    for (std::vector<uint>::const_iterator i = cupd->beginLReacs();
        i != cupd->endLReacs(); ++i)
    {
        s->computeReacProp(cidx, *i);
    }
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompMass(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));
    
    uint l_sidx = s->def()->comp(cidx)->specG2L(sidx);
    if (l_sidx == 0xFFFF) return 0.0;
    return (double)(s->fPoolCount[cidx][l_sidx]) / smath::AVOGADRO;
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompMass(State * s, uint cidx, uint sidx, double m)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));
    
    CompDef * cdef = s->def()->comp(cidx);
    uint l_sidx = cdef->specG2L(sidx);
    if (l_sidx == 0xFFFF) return;
    uint m2 = (uint)(m * smath::AVOGADRO);
    s->fPoolCount[cidx][l_sidx] = m2;
    
    // Recompute propensities.
    CompUpd * cupd = cdef->updateSpec(l_sidx);
    for (std::vector<uint>::const_iterator i = cupd->beginLReacs();
        i != cupd->endLReacs(); ++i)
    {
        s->computeReacProp(cidx, *i);
    }
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompConc(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));
    
    uint l_sidx = s->def()->comp(cidx)->specG2L(sidx);
    if (l_sidx == 0xFFFF) return 0.0;
    double vol = s->fCompVols[cidx];
    if (vol <= 0.0) return 0.0;
    return (double)(s->fPoolCount[cidx][l_sidx]) / 
        (1.0e3 * vol * smath::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompConc(State * s, uint cidx, uint sidx, double c)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx));
    assert(s->def()->isValidSpec(sidx));
    
    CompDef * cdef = s->def()->comp(cidx);
    uint l_sidx = cdef->specG2L(sidx);
    if (l_sidx == 0xFFFF) return;
    double vol = s->fCompVols[cidx];
    if (vol <= 0.0) return;
    uint m2 = (uint)(c * 1.0e3 * vol * smath::AVOGADRO);
    s->fPoolCount[cidx][l_sidx] = m2;
    
    // Recompute propensities.
    CompUpd * cupd = cdef->updateSpec(l_sidx);
    for (std::vector<uint>::const_iterator i = cupd->beginLReacs();
        i != cupd->endLReacs(); ++i)
    {
        s->computeReacProp(cidx, *i);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompClamped(State * s, uint cidx, uint sidx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidSpec(sidx) == true);
    
    uint l_sidx = s->def()->comp(cidx)->specG2L(sidx);
    if (l_sidx == 0xFFFF) return false;
    return ((s->fPoolFlags[cidx][l_sidx] & State::CLAMPED_POOLFLAG) != 0);
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompClamped(State * s, uint cidx, uint sidx, bool buf)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidSpec(sidx) == true);
    
    uint l_sidx = s->def()->comp(cidx)->specG2L(sidx);
    if (l_sidx == 0xFFFF) return;
    s->fPoolFlags[cidx][l_sidx] |= State::CLAMPED_POOLFLAG;
    if (buf == false) s->fPoolFlags[cidx][l_sidx] -= State::CLAMPED_POOLFLAG;
}

////////////////////////////////////////////////////////////////////////////////

double siGetCompReacK(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return 0.0;
    return s->fReacKcsts[cidx][l_ridx];
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacK(State * s, uint cidx, uint ridx, double kf)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return;
    assert(kf >= 0.0);
    
    s->fReacKcsts[cidx][l_ridx] = kf;
    s->computeCcst(cidx, l_ridx);
    s->computeReacProp(cidx, l_ridx);
}

////////////////////////////////////////////////////////////////////////////////

bool siGetCompReacActive(State * s, uint cidx, uint ridx)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return false;
    return ((s->fReacFlags[cidx][l_ridx] & State::ACTIVE_REACFLAG) != 0);
}

////////////////////////////////////////////////////////////////////////////////

void siSetCompReacActive(State * s, uint cidx, uint ridx, bool act)
{
    assert(s != 0);
    assert(s->def()->isValidComp(cidx) == true);
    assert(s->def()->isValidReac(ridx) == true);
    
    uint l_ridx = s->def()->comp(cidx)->reacG2L(ridx);
    if (l_ridx == 0xFFFF) return;
    s->fReacFlags[cidx][l_ridx] |= State::ACTIVE_REACFLAG;
    if (act == false) s->fReacFlags[cidx][l_ridx] -= State::ACTIVE_REACFLAG;
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

// END
