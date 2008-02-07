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

// Standard library headers.
#include <cassert>
#include <cmath>
#include <iostream>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/wmdirect/solver_core/comp.hpp>
#include <steps/wmdirect/solver_core/kproc.hpp>
#include <steps/wmdirect/solver_core/patch.hpp>
#include <steps/wmdirect/solver_core/sched.hpp>
#include <steps/wmdirect/solver_core/sreac.hpp>
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

SReac::SReac(ssim::SReacDef * srdef, Patch * patch)
: KProc()
, pSReacDef(srdef)
, pPatch(patch)
, pUpdVec()
, pCcst(0.0)
{
    assert(srdef != 0);
    assert(patch != 0);
    double vol;
    if (srdef->inside() == true)
    {
        assert(patch->iComp() != 0);
        vol = patch->iComp()->vol();
    }
    else
    {
        assert(patch->oComp() != 0);
        vol = patch->oComp()->vol();
    }
    pCcst = comp_ccst(srdef->kcst(), vol, srdef->order());
    assert(pCcst >= 0.0);
}
    
////////////////////////////////////////////////////////////////////////////////

SReac::~SReac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void SReac::setupDeps(void)
{
    Comp * icomp = pPatch->iComp();
    Comp * ocomp = pPatch->oComp();
    
    ssim::gidxTVecCI s_beg = def()->beginUpdColl_S();
    ssim::gidxTVecCI s_end = def()->endUpdColl_S();
    ssim::gidxTVecCI i_beg = def()->beginUpdColl_I();
    ssim::gidxTVecCI i_end = def()->endUpdColl_I();
    ssim::gidxTVecCI o_beg = def()->beginUpdColl_O();
    ssim::gidxTVecCI o_end = def()->endUpdColl_O();
    
    SchedIDXSet updset;
    
    KProcPVecCI kprocend = pPatch->kprocEnd();
    for (KProcPVecCI k = pPatch->kprocBegin(); k != kprocend; ++k)
    {
        for (ssim::gidxTVecCI spec = s_beg; spec != s_end; ++spec)
        {
            if ((*k)->depSpecPatch(*spec, pPatch) == true) 
                updset.insert((*k)->schedIDX());
        }
    }
    
    if (icomp != 0)
    {
        kprocend = icomp->kprocEnd();
        for (KProcPVecCI k = icomp->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
            {
                if ((*k)->depSpecComp(*spec, icomp) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
        
        PatchPVecCI ip_bgn = icomp->beginIPatches();
        PatchPVecCI ip_end = icomp->endIPatches();
        for (PatchPVecCI ip = ip_bgn; ip != ip_end; ++ip)
        {
            kprocend = (*ip)->kprocEnd();
            for (KProcPVecCI k = (*ip)->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
                {
                    if ((*k)->depSpecComp(*spec, icomp) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
        
        PatchPVecCI op_bgn = icomp->beginOPatches();
        PatchPVecCI op_end = icomp->endOPatches();
        for (PatchPVecCI op = op_bgn; op != op_end; ++op)
        {
            kprocend = (*op)->kprocEnd();
            for (KProcPVecCI k = (*op)->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
                {
                    if ((*k)->depSpecComp(*spec, icomp) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }
    
    if (ocomp != 0)
    {
        kprocend = ocomp->kprocEnd();
        for (KProcPVecCI k = ocomp->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
            {
                if ((*k)->depSpecComp(*spec, ocomp) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
        
        PatchPVecCI ip_bgn = ocomp->beginIPatches();
        PatchPVecCI ip_end = ocomp->endIPatches();
        for (PatchPVecCI ip = ip_bgn; ip != ip_end; ++ip)
        {
            kprocend = (*ip)->kprocEnd();
            for (KProcPVecCI k = (*ip)->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
                {
                    if ((*k)->depSpecComp(*spec, ocomp) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
        
        PatchPVecCI op_bgn = ocomp->beginOPatches();
        PatchPVecCI op_end = ocomp->endOPatches();
        for (PatchPVecCI op = op_bgn; op != op_end; ++op)
        {
            kprocend = (*op)->kprocEnd();
            for (KProcPVecCI k = (*op)->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
                {
                    if ((*k)->depSpecComp(*spec, ocomp) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }
    
    schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecComp(uint gidx, Comp * comp)
{
    if (comp == pPatch->iComp())
    {
        return (def()->dep_I(gidx) != ssim::DEP_NONE);
    }
    else if (comp == pPatch->oComp())
    {
        return (def()->dep_O(gidx) != ssim::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecPatch(uint gidx, Patch * patch)
{
    if (patch != pPatch) return false;
    return (def()->dep_S(gidx) != ssim::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

void SReac::reset(void)
{
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

double SReac::rate(void) const
{    
    if (inactive()) return 0.0;
    
    // First we compute the combinatorial part.
    //   1/ for the surface part of the stoichiometry
    //   2/ for the inner or outer volume part of the stoichiometry,
    //      depending on whether the sreac is inner() or outer()
    // Then we multiply with mesoscopic constant.
    
    ssim::PatchDef * pdef = pPatch->def();
    ssim::lidxT lidx = pdef->sreacG2L(def()->gidx());
    
    double h_mu = 1.0;
    
    uint * lhs_s_vec = pdef->sreac_lhs_S_bgn(lidx);
    uint * cnt_s_vec = pPatch->pools();
    uint nspecs_s = pdef->countSpecs();
    for (uint s = 0; s < nspecs_s; ++s)
    {
        uint lhs = lhs_s_vec[s];
        if (lhs == 0) continue;
        uint cnt = cnt_s_vec[s];
        if (lhs > cnt) 
        {
            return 0.0;
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
    
    if (def()->inside())
    {
        uint * lhs_i_vec = pdef->sreac_lhs_I_bgn(lidx);
        uint * cnt_i_vec = pPatch->iComp()->pools();
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            uint lhs = lhs_i_vec[s];
            if (lhs == 0) continue;
            uint cnt = cnt_i_vec[s];
            if (lhs > cnt) 
            {
                return 0.0;
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
    }
    else if (def()->outside())
    {
        uint * lhs_o_vec = pdef->sreac_lhs_O_bgn(lidx);
        uint * cnt_o_vec = pPatch->oComp()->pools();
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            uint lhs = lhs_o_vec[s];
            if (lhs == 0) continue;
            uint cnt = cnt_o_vec[s];
            if (lhs > cnt) 
            {
                return 0.0;
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
    }
    
    return h_mu * pCcst;
}

////////////////////////////////////////////////////////////////////////////////

SchedIDXVec const & SReac::apply(State * state)
{   
    ssim::PatchDef * pdef = pPatch->def();
    ssim::lidxT lidx = pdef->sreacG2L(def()->gidx());
    
    // Update patch pools.
    int * upd_s_vec = pdef->sreac_upd_S_bgn(lidx);
    uint * cnt_s_vec = pPatch->pools();
    uint nspecs_s = pdef->countSpecs();
    
    for (uint s = 0; s < nspecs_s; ++s)
    {
        if (pPatch->clamped(s) == true) continue;
        int upd = upd_s_vec[s];
        if (upd == 0) continue;
        int nc = static_cast<int>(cnt_s_vec[s]) + upd; 
        assert(nc >= 0);
        cnt_s_vec[s] = static_cast<uint>(nc);
    }
    
    // Update inner comp pools.
    Comp * icomp = pPatch->iComp();
    if (icomp != 0)
    {
        int * upd_i_vec = pdef->sreac_upd_I_bgn(lidx);
        uint * cnt_i_vec = icomp->pools();
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            if (icomp->clamped(s) == true) continue;
            int upd = upd_i_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            assert(nc >= 0);
            cnt_i_vec[s] = static_cast<uint>(nc);
        }
    }
    
    // Update outer comp pools.
    Comp * ocomp = pPatch->oComp();
    if (ocomp != 0)
    {
        int * upd_o_vec = pdef->sreac_upd_O_bgn(lidx);
        uint * cnt_o_vec = ocomp->pools();
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            if (ocomp->clamped(s) == true) continue;
            int upd = upd_o_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            assert(nc >= 0);
            cnt_o_vec[s] = static_cast<uint>(nc);
        }
    }
    
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
