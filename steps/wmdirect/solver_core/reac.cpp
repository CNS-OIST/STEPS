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

// Standard library & STL headers.
#include <cmath>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/wmdirect/solver_core/comp.hpp>
#include <steps/wmdirect/solver_core/kproc.hpp>
#include <steps/wmdirect/solver_core/patch.hpp>
#include <steps/wmdirect/solver_core/reac.hpp>
#include <steps/wmdirect/solver_core/sched.hpp>
#include <steps/wmdirect/solver_core/state.hpp>

NAMESPACE_ALIAS(steps::math, smath);
NAMESPACE_ALIAS(steps::sim, ssim);

////////////////////////////////////////////////////////////////////////////////

static inline double comp_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(ssim::ReacDef * rdef, Comp * comp)
: KProc()
, pReacDef(rdef)
, pComp(comp)
, pUpdVec()
, pCcst()
{
    assert(pReacDef != 0);
    assert(pComp != 0);
    pCcst = comp_ccst(pReacDef->kcst(), pComp->vol(), pReacDef->order());
    assert(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setupDeps(void)
{
    SchedIDXSet updset; 
    ssim::gidxTVecCI sbgn = def()->bgnUpdColl();
    ssim::gidxTVecCI send = def()->endUpdColl();
    
    // Search in local compartment.
    KProcPVecCI kprocend = pComp->kprocEnd();
    for (KProcPVecCI k = pComp->kprocBegin(); k != kprocend; ++k)
    {
        for (ssim::gidxTVecCI s = sbgn; s != send; ++s)
        {
            if ((*k)->depSpecComp(*s, pComp) == true) 
                updset.insert((*k)->schedIDX());
        }
    }
    
    // Search in neighbouring patches.
    for (PatchPVecCI p = pComp->beginIPatches(); p != pComp->endIPatches(); ++p)
    {        
        kprocend = (*p)->kprocEnd();
        for (KProcPVecCI k = (*p)->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecComp(*s, pComp) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
    }
    for (PatchPVecCI p = pComp->beginOPatches(); p != pComp->endOPatches(); ++p)
    {        
        kprocend = (*p)->kprocEnd();
        for (KProcPVecCI k = (*p)->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI s = sbgn; s != send; ++s)
            {
                if ((*k)->depSpecComp(*s, pComp) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
    }
    
    schedIDXSet_To_Vec(updset, pUpdVec);   
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecComp(uint gidx, Comp * comp)
{
    if (pComp != comp) return false;
    return def()->dep(gidx);
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecPatch(uint gidx, Patch * tri)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void Reac::reset(void)
{
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

double Reac::rate(void) const
{
    if (inactive()) return 0.0;
    
    // Prefetch some variables.
    ssim::CompDef * cdef = pComp->def();
    uint nspecs = cdef->countSpecs();
    uint * lhs_vec = cdef->reac_lhs_bgn(cdef->reacG2L(def()->gidx()));
    uint * cnt_vec = pComp->pools();
    
    // Compute combinatorial part.
    double h_mu = 1.0;
    for (uint pool = 0; pool < nspecs; ++pool)
    {
        uint lhs = lhs_vec[pool];
        if (lhs == 0) continue;
        uint cnt = cnt_vec[pool];
        if (lhs > cnt) 
        {
            h_mu = 0.0;
            break;
        }
        switch (lhs)
        {
            case 4:
            {
                h_mu *= static_cast<double>(cnt - 3);
            }
            case 3:
            {
                h_mu *= static_cast<double>(cnt - 2);
            }
            case 2:
            {
                h_mu *= static_cast<double>(cnt - 1);
            }
            case 1:
            {
                h_mu *= static_cast<double>(cnt);
                break;
            }
            default:
            {
                assert(0);
                return 0.0;
            }
        }
    }

    // Multiply with scaled reaction constant.
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

SchedIDXVec const & Reac::apply(State * s)
{
    uint * local = pComp->pools();
    ssim::CompDef * cdef = pComp->def();
    uint l_ridx = cdef->reacG2L(def()->gidx());
    int * upd_vec = cdef->reac_upd_bgn(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        if (pComp->clamped(i) == true) continue;
        int j = upd_vec[i];
        if (j == 0) continue;
        int nc = static_cast<int>(local[i]) + j;
        local[i] = static_cast<uint>(nc);
    }
    rExtent++;
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
