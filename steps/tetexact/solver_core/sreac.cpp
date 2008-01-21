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
#include <cassert>
#include <cmath>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/types.hpp>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>
#include <steps/tetexact/solver_core/sreac.hpp>
#include <steps/tetexact/solver_core/tet.hpp>
#include <steps/tetexact/solver_core/tri.hpp>

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

SReac::SReac(ssim::SReacDef * srdef, Tri * tri)
: KProc()
, pSReacDef(srdef)
, pTri(tri)
, pUpdVec()
, pCcst(0.0)
{
    assert(srdef != 0);
    assert(tri != 0);
    double vol;
    if (srdef->inside() == true)
    {
        assert(tri->iTet() != 0);
        vol = tri->iTet()->vol();
    }
    else
    {
        assert(tri->oTet() != 0);
        vol = tri->oTet()->vol();
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
    // For all non-zero entries gidx in SReacDef's UPD_S:
    //   Perform depSpecTri(gidx,tri()) for:
    //     All kproc's of tri()
    //     NOTE: we currently don't test kproc's of inner and external tet, 
    //     because that's not strictly necessary.
    //
    // If inner tetrahedron exists:
    //   For all non-zero entries gidx in SReacDef's UPD_I:
    //     Perform depSpecTet(gidx,itet) for:
    //       All kproc's of itet
    //       All kproc's of triangles next to itet
    //
    // If outer tetrahedron exists:
    //   Similar to inner tet.
    //
    // All dependencies are first collected into a std::set, to sort them
    // and to eliminate duplicates. At the end of the routine, they are
    // copied into the vector that will be returned during execution.
    
    Tet * itet = tri()->iTet();
    Tet * otet = tri()->oTet();
    
    ssim::gidxTVecCI s_beg = def()->beginUpdColl_S();
    ssim::gidxTVecCI s_end = def()->endUpdColl_S();
    ssim::gidxTVecCI i_beg = def()->beginUpdColl_I();
    ssim::gidxTVecCI i_end = def()->endUpdColl_I();
    ssim::gidxTVecCI o_beg = def()->beginUpdColl_O();
    ssim::gidxTVecCI o_end = def()->endUpdColl_O();
    
    SchedIDXSet updset;
    
    KProcPVecCI kprocend = tri()->kprocEnd();
    for (KProcPVecCI k = tri()->kprocBegin(); k != kprocend; ++k)
    {
        for (ssim::gidxTVecCI spec = s_beg; spec != s_end; ++spec)
        {
            if ((*k)->depSpecTri(*spec, tri()) == true) 
                updset.insert((*k)->schedIDX());
        }
    }
    
    if (itet != 0)
    {
        kprocend = itet->kprocEnd();
        for (KProcPVecCI k = itet->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
            {
                if ((*k)->depSpecTet(*spec, itet) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
        
        for (uint i = 0; i < 4; ++i)
        {
            Tri * tri = itet->nextTri(i);
            if (tri == 0) continue;
            kprocend = tri->kprocEnd();
            for (KProcPVecCI k = tri->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = i_beg; spec != i_end; ++spec)
                {
                    if ((*k)->depSpecTet(*spec, itet) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }
    
    if (otet != 0)
    {
        kprocend = otet->kprocEnd();
        for (KProcPVecCI k = otet->kprocBegin(); k != kprocend; ++k)
        {
            for (ssim::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
            {
                if ((*k)->depSpecTet(*spec, otet) == true) 
                    updset.insert((*k)->schedIDX());
            }
        }
        
        for (uint i = 0; i < 4; ++i)
        {
            Tri * tri = otet->nextTri(i);
            if (tri == 0) continue;
            kprocend = tri->kprocEnd();
            for (KProcPVecCI k = tri->kprocBegin(); k != kprocend; ++k)
            {
                for (ssim::gidxTVecCI spec = o_beg; spec != o_end; ++spec)
                {
                    if ((*k)->depSpecTet(*spec, otet) == true) 
                        updset.insert((*k)->schedIDX());
                }
            }
        }
    }
    
    schedIDXSet_To_Vec(updset, pUpdVec);
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecTet(uint gidx, Tet * tet)
{
    // We need to check whether the tet is inside or outside.
    //   -> If inside: check dependency using SReacDef's I_DEP
    //   -> If outside: check dependency using SReacDef's O_DEP
    //   -> If neither, return.
    //
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp
    
    if (tet == tri()->iTet())
    {
        return (def()->dep_I(gidx) != ssim::DEP_NONE);
    }
    else if (tet == tri()->oTet())
    {
        return (def()->dep_O(gidx) != ssim::DEP_NONE);
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReac::depSpecTri(uint gidx, Tri * triangle)
{
    if (triangle != tri()) return false;
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
    
    ssim::PatchDef * pdef = tri()->patchdef();
    ssim::lidxT lidx = pdef->sreacG2L(def()->gidx());
    
    double h_mu = 1.0;
    
    uint * lhs_s_vec = pdef->sreac_lhs_S_bgn(lidx);
    uint * cnt_s_vec = tri()->pools();
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
        uint * cnt_i_vec = tri()->iTet()->pools();
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
        uint * cnt_o_vec = tri()->oTet()->pools();
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

SchedIDXVec const & SReac::apply(State * s)
{
    ssim::PatchDef * pdef = tri()->patchdef();
    ssim::lidxT lidx = pdef->sreacG2L(def()->gidx());
    
    // Update triangle pools.
    int * upd_s_vec = pdef->sreac_upd_S_bgn(lidx);
    uint * cnt_s_vec = tri()->pools();
    uint nspecs_s = pdef->countSpecs();
    for (uint s = 0; s < nspecs_s; ++s)
    {
        int upd = upd_s_vec[s];
        if (upd == 0) continue;
        int nc = static_cast<int>(cnt_s_vec[s]) + upd;
        assert(nc > 0);
        cnt_s_vec[s] = static_cast<uint>(nc);
    }
    
    // Update inner tet pools.
    Tet * itet = tri()->iTet();
    if (itet != 0)
    {
        int * upd_i_vec = pdef->sreac_upd_I_bgn(lidx);
        uint * cnt_i_vec = itet->pools();
        uint nspecs_i = pdef->countSpecs_I();
        for (uint s = 0; s < nspecs_i; ++s)
        {
            int upd = upd_i_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_i_vec[s]) + upd;
            assert(nc > 0);
            cnt_i_vec[s] = static_cast<uint>(nc);
        }
    }
    
    // Update outer tet pools.
    Tet * otet = tri()->oTet();
    if (otet != 0)
    {
        int * upd_o_vec = pdef->sreac_upd_O_bgn(lidx);
        uint * cnt_o_vec = otet->pools();
        uint nspecs_o = pdef->countSpecs_O();
        for (uint s = 0; s < nspecs_o; ++s)
        {
            int upd = upd_o_vec[s];
            if (upd == 0) continue;
            int nc = static_cast<int>(cnt_o_vec[s]) + upd;
            assert(nc > 0);
            cnt_o_vec[s] = static_cast<uint>(nc);
        }
    }
    
    // Update outer tet pools.
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
