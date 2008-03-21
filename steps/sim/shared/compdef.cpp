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
#include <iostream>
#include <set>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>

USING(std, string);
USING(std, vector);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

CompDef::CompDef(StateDef * sdef, gidxT idx, string const & name)
: pStateDef(sdef)
, pGIDX(idx)
, pName(name)
, pVolume(0.0)
, pLocalIndicesSetupDone(false)
, pIPatches()
, pOPatches()
, pSpec_N(0)
, pSpec_G2L(0)
, pSpec_L2G(0)
, pReac_N(0)
, pReac_G2L(0)
, pReac_L2G(0)
, pReac_DEP_Spec(0)
, pReac_LHS_Spec(0)
, pReac_UPD_Spec(0)
, pDiff_N(0)
, pDiff_G2L(0)
, pDiff_L2G(0)
, pDiff_DEP_Spec(0)
, pDiff_LIG(0)
{
    uint nspecs = statedef()->countSpecs();
    pSpec_G2L = new lidxT[nspecs];
    std::fill_n(pSpec_G2L, nspecs, LIDX_UNDEFINED);
    uint nreacs = statedef()->countReacs();
    pReac_G2L = new lidxT[nreacs];
    std::fill_n(pReac_G2L, nreacs, LIDX_UNDEFINED);
    uint ndiffs = statedef()->countDiffs();
    pDiff_G2L = new lidxT[ndiffs];
    std::fill_n(pDiff_G2L, ndiffs, LIDX_UNDEFINED);
}

////////////////////////////////////////////////////////////////////////////////

CompDef::~CompDef(void)
{
    delete[] pSpec_G2L;
    delete[] pSpec_L2G;
    delete[] pReac_G2L;
    delete[] pReac_L2G;
    delete[] pReac_DEP_Spec;
    delete[] pReac_LHS_Spec;
    delete[] pReac_UPD_Spec;
    delete[] pDiff_G2L;
    delete[] pDiff_L2G;
    delete[] pDiff_DEP_Spec;
    delete[] pDiff_LIG;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addSpec(gidxT idx)
{
    assert(pLocalIndicesSetupDone == false);
    assert(statedef()->spec(idx) != 0);
    if (pSpec_G2L[idx] != LIDX_UNDEFINED) return;
    pSpec_G2L[idx] = pSpec_N++;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addReac(gidxT idx)
{
    assert(pLocalIndicesSetupDone == false);
    assert(statedef()->reac(idx) != 0);
    if (pReac_G2L[idx] != LIDX_UNDEFINED) return;
    pReac_G2L[idx] = pReac_N++;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addDiff(gidxT idx)
{
    assert(pLocalIndicesSetupDone == false);
    assert(statedef()->diff(idx) != 0);
    if (pDiff_G2L[idx] != LIDX_UNDEFINED) return;
    pDiff_G2L[idx] = pDiff_N++;
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addReferences(void)
{
    uint ngspecs = statedef()->countSpecs();
    uint ngreacs = statedef()->countReacs();
    uint ngdiffs = statedef()->countDiffs();
    
    for (uint r = 0; r < ngreacs; ++r)
    {
        if (pReac_G2L[r] == LIDX_UNDEFINED) continue;
        ReacDef * rdef = statedef()->reac(r);
        assert(rdef != 0);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (rdef->req(s) == true) addSpec(s);
        }
    }
    
    for (uint d = 0; d < ngdiffs; ++d)
    {
        if (pDiff_G2L[d] == LIDX_UNDEFINED) continue;
        DiffDef * ddef = statedef()->diff(d);
        assert(ddef != 0);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (ddef->req(s) == true) addSpec(s);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::setupLocalIndices(void)
{
    assert(pLocalIndicesSetupDone == false);
    pLocalIndicesSetupDone = true;
    
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
    
    uint ngspecs = statedef()->countSpecs();
    
    if (pReac_N != 0)
    {
        pReac_L2G = new gidxT[pReac_N];
        uint ngreacs = statedef()->countReacs();
        for (uint i = 0; i < ngreacs; ++i)
        {
            lidxT lidx = reacG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pReac_L2G[lidx] = i;
        }
        
        uint arrsize = countSpecs() * countReacs();
        pReac_DEP_Spec = new depT[arrsize];
        pReac_LHS_Spec = new uint[arrsize];
        pReac_UPD_Spec = new int[arrsize];
        // DEBUG: 08-Feb-2008
        std::fill_n(pReac_DEP_Spec, arrsize, 0);
        std::fill_n(pReac_LHS_Spec, arrsize, 0);
        std::fill_n(pReac_UPD_Spec, arrsize, 0);
        for (uint ri = 0; ri < countReacs(); ++ri)
        {
            ReacDef * reacdef = reac(ri);
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (reacdef->req(si) == false) continue;
                
                // TODO: turn into error check?
                lidxT sil = specG2L(si);
                assert(sil != LIDX_UNDEFINED);
                
                uint aridx = _IDX_Reac_Spec(ri, sil);
                pReac_DEP_Spec[aridx] = reacdef->dep(si);
                pReac_LHS_Spec[aridx] = reacdef->lhs(si);
                pReac_UPD_Spec[aridx] = reacdef->upd(si);
            }
        }
    }
    
    if (pDiff_N != 0)
    {
        pDiff_L2G = new gidxT[pDiff_N];
        uint ngdiffs = statedef()->countDiffs();
        for (uint i = 0; i < ngdiffs; ++i)
        {
            lidxT lidx = diffG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pDiff_L2G[lidx] = i;
        }
        
        uint arrsize = countSpecs() * countDiffs();
        pDiff_DEP_Spec = new depT[arrsize];
        // DEBUG: 08-Feb-2008
        std::fill_n(pReac_DEP_Spec, arrsize, 0);
        pDiff_LIG = new lidxT[countDiffs()];
        for (uint di = 0; di < countDiffs(); ++di)
        {
            DiffDef * diffdef = diff(di);
            pDiff_LIG[di] = specG2L(diffdef->lig());
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (diffdef->req(si) == false) continue;
                
                // TODO: turn into error check?
                lidxT sil = specG2L(si);
                assert(sil != LIDX_UNDEFINED);
                
                uint aridx = _IDX_Diff_Spec(di, sil);
                pDiff_DEP_Spec[aridx] = diffdef->dep(di);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::setupDependencies(void)
{
    /*
    // Prefetch a couple of useful variables.
    uint nspecs = countSpecs();
    uint nreacs = countReacs();
    uint ndiffs = countDiffs();
    
    // Create a list of CompUpd objects, one for each species variable,
    // that summarize which (local) processes and rules need to be 
    // updated when the species is changed.
    pSpecUpd = new CompUpd[nspecs];
    
    // In the following loop, a number of things are set up:
    //   - for each reaction, its stoichiometric dependencies and update 
    //     vectors. These are stored as C arrays of signed and unsigned
    //     integers, respectively.
    //   - add the reaction to the CompUpd objects of all variables on
    //     which the reaction is dependent.
    pReacSpecDeps = new uint[nspecs * nreacs];
    pReacSpecUpds = new int[nspecs * nreacs];
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
    {
        ReacDef * rdef = reac(l_ridx);
        uint idx = l_ridx * nspecs;
        for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
        {
            uint g_sidx = specL2G(l_sidx);
            uint lhs = rdef->lhs(g_sidx);
            uint rhs = rdef->rhs(g_sidx);
            pReacSpecDeps[idx] = lhs;
            pReacSpecUpds[idx] = static_cast<int>(rhs) - static_cast<int>(lhs);
            if (lhs > 0) pSpecUpd[l_sidx].pLReacs.push_back(l_ridx);
            idx++;
        }
    }
    
    // For each diffusion rule, add it to the updates for the ligand
    // species.
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
    {
        DiffDef * ddef = diff(l_didx);
        uint l_sidx = diffG2L(ddef->lig());
        pSpecUpd[l_sidx].pLDiffs.push_back(l_didx);
    }
    
    // For each reaction, build the set of required updates by merging
    // the CompUpd objects of each species that is non-zero in the
    // reaction's stoichiometric update vector.
    pReacUpd = new CompUpd[nreacs];
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
    {
        uint idx = l_ridx * nspecs;
        for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
        {
            if (pReacSpecUpds[idx] != 0)
                pReacUpd[l_ridx].merge(pSpecUpd[l_sidx]);
            idx++;
        }
    }
    
    // For each diffusion, build the set of required updates by copying
    // the CompUpd object of the ligand.
    pDiffUpd = new CompUpd[ndiffs];
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
    {
        DiffDef * ddef = diff(l_didx);
        uint l_sidx = diffG2L(ddef->lig());
        pDiffUpd[l_didx].merge(pSpecUpd[l_sidx]);
    }
    
    // Perform compacting step on all CompUpd objects (just to be sure).
    for (uint l_sidx = 0; l_sidx < nspecs; ++l_sidx)
        pSpecUpd[l_sidx].compact();
    for (uint l_ridx = 0; l_ridx < nreacs; ++l_ridx)
        pSpecUpd[l_ridx].compact();
    for (uint l_didx = 0; l_didx < ndiffs; ++l_didx)
        pSpecUpd[l_didx].compact();
    */
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addIPatchDef(PatchDef * p)
{
    // Make some checks.
    assert(p != 0);
    assert(p->ocompdef() == this);
    // Check whether it's already included.
    PatchDefPVecI ip_end = pIPatches.end();
    if (std::find(pIPatches.begin(), ip_end, p) != ip_end) return;
#ifndef NDEBUG
    PatchDefPVecI op_end = pOPatches.end();
    assert(std::find(pOPatches.begin(), op_end, p) == op_end);
#endif
    // Include.
    pIPatches.push_back(p);
}

////////////////////////////////////////////////////////////////////////////////

void CompDef::addOPatchDef(PatchDef * p)
{
    // Make some checks.
    assert(p != 0);
    assert(p->icompdef() == this);
    // Check whether it's already included.
    PatchDefPVecI op_end = pOPatches.end();
    if (std::find(pOPatches.begin(), op_end, p) != op_end) return;
#ifndef NDEBUG
    PatchDefPVecI ip_end = pIPatches.end();
    assert(std::find(pIPatches.begin(), ip_end, p) == ip_end);
#endif
    // Include.
    pOPatches.push_back(p);
}
    
////////////////////////////////////////////////////////////////////////////////

SpecDef * CompDef::spec(lidxT idx) const
{
    assert(pLocalIndicesSetupDone == true);
    return statedef()->spec(specL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

depT CompDef::reac_dep(lidxT reac, lidxT spec) const
{
    return pReac_DEP_Spec[spec + (reac * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

uint * CompDef::reac_lhs_bgn(lidxT reac) const
{
    return pReac_LHS_Spec + (reac * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

uint * CompDef::reac_lhs_end(lidxT reac) const
{
    return pReac_LHS_Spec + ((reac + 1) * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int * CompDef::reac_upd_bgn(lidxT reac) const
{
    return pReac_UPD_Spec + (reac * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int * CompDef::reac_upd_end(lidxT reac) const
{
    return pReac_UPD_Spec + ((reac + 1) * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

ReacDef * CompDef::reac(lidxT idx) const
{
    assert(pLocalIndicesSetupDone == true);
    return statedef()->reac(reacL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

depT CompDef::diff_dep(lidxT diff, lidxT spec) const
{
    return pDiff_DEP_Spec[spec + (diff * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

lidxT CompDef::diff_lig(lidxT diff) const
{
    return pDiff_LIG[diff];
}

////////////////////////////////////////////////////////////////////////////////

DiffDef * CompDef::diff(lidxT idx) const
{
    assert(pLocalIndicesSetupDone == true);
    return statedef()->diff(diffL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

// END
