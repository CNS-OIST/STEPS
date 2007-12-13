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

// STL headers.
#include <algorithm>
#include <cassert>
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>

USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

PatchDef::PatchDef
(
    StateDef * sdef, gidxT idx, 
    std::string const & name, 
    CompDef * inner, CompDef * outer
)
: pStateDef(sdef)
, pGIDX(idx)
, pName(name)
, pArea(0.0)
, pInner(inner)
, pOuter(outer)
, pSpec_N(0)
, pSpec_G2L(0)
, pSpec_L2G(0)
, pSReac_N(0)
, pSReac_G2L(0)
, pSReac_L2G(0)
, pSReac_DEP_I_Spec(0)
, pSReac_DEP_S_Spec(0)
, pSReac_DEP_O_Spec(0)
, pSReac_LHS_I_Spec(0)
, pSReac_LHS_S_Spec(0)
, pSReac_LHS_O_Spec(0)
, pSReac_UPD_I_Spec(0)
, pSReac_UPD_S_Spec(0)
, pSReac_UPD_O_Spec(0)
{
    assert(pStateDef != 0);
    if (pInner != 0)
    {
        assert(pStateDef == pInner->statedef());
    }
    if (pOuter != 0)
    {
        assert(pStateDef == pOuter->statedef());
    }
    
    // Pre-generate certain arrays.
    uint nspecs = pStateDef->countSpecs();
    pSpec_G2L = new lidxT[nspecs];
    std::fill_n(pSpec_G2L, nspecs, LIDX_UNDEFINED);
    uint nsreacs = pStateDef->countSReacs(); 
    pSReac_G2L = new lidxT[nsreacs]; 
    std::fill_n(pSReac_G2L, nsreacs, LIDX_UNDEFINED);
}

////////////////////////////////////////////////////////////////////////////////

PatchDef::~PatchDef(void)
{
    delete[] pSpec_G2L;
    delete[] pSpec_L2G;
    delete[] pSReac_G2L;
    delete[] pSReac_L2G;
    delete[] pSReac_DEP_I_Spec;
    delete[] pSReac_DEP_S_Spec;
    delete[] pSReac_DEP_O_Spec;
    delete[] pSReac_LHS_I_Spec;
    delete[] pSReac_LHS_S_Spec;
    delete[] pSReac_LHS_O_Spec;
    delete[] pSReac_UPD_I_Spec;
    delete[] pSReac_UPD_S_Spec;
    delete[] pSReac_UPD_O_Spec;
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::addSpec(gidxT idx)
{
    assert(statedef()->spec(idx) != 0);
    if (pSpec_G2L[idx] != LIDX_UNDEFINED) return;
    pSpec_G2L[idx] = pSpec_N++;
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::addSReac(gidxT idx)
{
    assert(statedef()->sreac(idx) != 0);
    if (pSReac_G2L[idx] != LIDX_UNDEFINED) return;
    pSReac_G2L[idx] = pSReac_N++;
}

////////////////////////////////////////////////////////////////////////////////
    
void PatchDef::setupLocalIndices(void)
{
    if (pSpec_N != 0)
    {
        pSpec_L2G = new gidxT[pSpec_N];
        uint ngspecs = statedef()->countSpecs();
        for (uint i = 0; i < ngspecs; ++i)
        {
            lidxT lidx = specG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pSpec_L2G[lidx] = i;
        }
    }
    
    if (pSReac_N != 0)
    {
        pSReac_L2G = new gidxT[pSReac_N];
        uint ngsreacs = statedef()->countSReacs();
        for (uint i = 0; i < ngsreacs; ++i)
        {
            lidxT lidx = sreacG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pSReac_L2G[lidx] = i;
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::setupDependencies(void)
{
    
}

////////////////////////////////////////////////////////////////////////////////

SpecDef * PatchDef::spec(lidxT idx) const
{
    return statedef()->spec(specL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

SReacDef * PatchDef::sreac(lidxT idx) const
{
    return statedef()->sreac(sreacL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

// END
