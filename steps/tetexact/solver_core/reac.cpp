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
#include <steps/sim/shared/compdef.hpp> 
#include <steps/sim/shared/reacdef.hpp>
#include <steps/tetexact/solver_core/reac.hpp>
#include <steps/tetexact/solver_core/state.hpp>
#include <steps/tetexact/solver_core/tet.hpp>

NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

inline double comp_ccst(double kcst, double vol, uint order)
{
    double vscale = 1.0e3 * vol * smath::AVOGADRO;
    int o1 = static_cast<int>(order) - 1;
    if (o1 < 0) o1 = 0;
    return kcst * pow(vscale, static_cast<double>(-o1));
}

////////////////////////////////////////////////////////////////////////////////

Reac::Reac(ReacDef * rdef, Tet * tet)
: pReacDef(rdef)
, pTet(tet)
, pUpdVec()
, pCcst(0.0)
{
    assert(pReacDef != 0);
    assert(pTet != 0);
    pCcst = comp_ccst(pReacDef->kcst(), pTet->vol(), pReacDef->order());
    assert(pCcst >= 0);
}

////////////////////////////////////////////////////////////////////////////////

Reac::~Reac(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void Reac::setupDeps(void)
{
    typedef std::vector<uint> uintVec;
    typedef uintVec::const_iterator uintVecCI;
    
    // Collect the global indices of all species that are updated by a 
    // single occurence of this reaction.
    CompDef * cdef = pTet->compdef();
    uint l_ridx = cdef->reacG2L(def()->gidx());
    int * lupds = cdef->reacSpecUpds(l_ridx);
    uint nspecs = cdef->countSpecs();
    uintVec upd;
    for (uint i = 0; i < nspecs; ++i)
    {
        if (lupds[i] != 0) upd.push_back(cdef->specL2G(i));
    }
    uintVecCI upd_end = upd.end();
    
    // Search for dependencies among kinetic processes in local voxel.
    KProcPVecCI kprocend = pTet->kprocEnd();
    for (KProcPVecCI k = pTet->kprocBegin(); k != kprocend; ++k)
    {
        for (uintVecCI u = upd.begin(); u != upd_end; ++u)
        {
            if ((*k)->depSpecTet(*u, pTet) == true) 
                pUpdVec.push_back((*k)->schedIDX());
        }
    }
    
    // Search for dependencies among kprocs in neighbouring tetrahedrons.
    // Strictly, this would currently never have to be the case...?
    for (uint i = 0; i < 4; ++i)
    {
        // Fetch next tetrahedron, if it exists.
        Tet * next = pTet->nextTet(i);
        if (next == 0) continue;
        
        kprocend = next->kprocEnd();
        for (KProcPVecCI k = next->kprocBegin(); k != kprocend; ++k)
        {
            for (uintVecCI u = upd.begin(); u != upd_end; ++u)
            {
                if ((*k)->depSpecTet(*u, pTet) == true) 
                    pUpdVec.push_back((*k)->schedIDX());
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Reac::depSpecTet(uint gidx, Tet * tet)
{
    if (pTet != tet) return false;
    return def()->dependsOnSpec(gidx);
}

////////////////////////////////////////////////////////////////////////////////

void Reac::reset(void)
{
}

////////////////////////////////////////////////////////////////////////////////

double Reac::rate(void) const
{
    // Prefetch some variables.
    CompDef * cdef = pTet->compdef();
    uint nspecs = cdef->countSpecs();
    uint * deps = cdef->reacSpecDeps(cdef->reacG2L(def()->gidx()));
    uint * cnts = pTet->pools();
    
    // Compute combinatorial part.
    double h_mu = 1.0;
    for (uint pool = 0; pool < nspecs; ++pool)
    {
        uint dep = deps[pool];
        if (dep == 0) continue;
        uint cnt = cnts[pool];
        if (dep > cnt) 
        {
            h_mu = 0.0;
            break;
        }
        switch (dep)
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
    uint * local = pTet->pools();
    CompDef * cdef = pTet->compdef();
    uint l_ridx = cdef->reacG2L(def()->gidx());
    int * upds = cdef->reacSpecUpds(l_ridx);
    uint nspecs = cdef->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        int j = upds[i];
        if (j == 0) continue;
        local[i] += j;
    }
    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
