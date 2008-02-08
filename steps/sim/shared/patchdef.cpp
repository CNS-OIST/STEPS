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
, pLocalIndicesSetupDone(false)
, pSpec_N_I(0)
, pSpec_N_S(0)
, pSpec_N_O(0)
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
        pInner->addOPatchDef(this);
    }
    if (pOuter != 0)
    {
        assert(pStateDef == pOuter->statedef());
        pOuter->addIPatchDef(this);
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
    assert(pLocalIndicesSetupDone == false);
    assert(statedef()->spec(idx) != 0);
    if (pSpec_G2L[idx] != LIDX_UNDEFINED) return;
    pSpec_G2L[idx] = pSpec_N_S++;
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::addSReac(gidxT idx)
{
    assert(pLocalIndicesSetupDone == false);
    assert(statedef()->sreac(idx) != 0);
    if (pSReac_G2L[idx] != LIDX_UNDEFINED) return;
    pSReac_G2L[idx] = pSReac_N++;
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::addReferences(void)
{
    uint ngspecs = statedef()->countSpecs();
    uint ngsreacs = statedef()->countSReacs();
    
    for (uint sr = 0; sr < ngsreacs; ++sr)
    {
        if (pSReac_G2L[sr] == LIDX_UNDEFINED) continue;
        SReacDef * srdef = statedef()->sreac(sr);
        for (uint s = 0; s < ngspecs; ++s)
        {
            if (srdef->req_S(s) == true) addSpec(s);
            if (srdef->req_I(s) == true)
            {
                assert(icompdef() != 0);
                icompdef()->addSpec(s);
            }
            if (srdef->req_O(s) == true)
            {
                assert(ocompdef() != 0);
                ocompdef()->addSpec(s);
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////
    
void PatchDef::setupLocalIndices(void)
{
    // 1 -- DEAL WITH PATCH SPECIES
    // (Only if any species have been added to the patch)
    //   -> Setup local indices for all species.
    //
    // 2 -- COPY #SPECIES FOR INNER AND OUTER COMPS
    //   (These are required a lot during simulation, so it's sound
    //   to have them ready here to avoid an extra level of pointer 
    //   lookup.)
    //
    // 3 -- DEAL WITH PATCH SREAC'S
    // (Only if any surface reactions have been added to the patch)
    //   -> Setup local indices for all surface reactions.
    //   -> The SReac objects have LHS, DEP and UPD vectors expressed in
    //      global species indices. Transform this to local indices:
    //      -> Pre-create pSReac_DEP, _LHS and _UPD vectors
    //          -> Always for surface (_S)
    //          -> For inner comp (_I) if inner comp has been defined
    //          -> For outer comp (_O) if outer comp has been defined 
    //      -> Loop over the SReacDef objects added to this patch:
    //          -> Fill out the newly created vectors by appropriate
    //             copying of the vectors defined the SReacDef object.
    //          -> While doing this, check whether everything can be 
    //             resolved.
    
    assert(pLocalIndicesSetupDone == false);
    pLocalIndicesSetupDone = true;
        
    // 1 -- DEAL WITH PATCH SPECIES
    uint ngspecs = statedef()->countSpecs();
    if (countSpecs() != 0) 
    {
        pSpec_L2G = new gidxT[countSpecs()];
        for (uint i = 0; i < ngspecs; ++i)
        {
            lidxT lidx = specG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pSpec_L2G[lidx] = i;
        }
    }
    
    // 2 -- COPY #SPECS FOR INNER AND OUTER COMPS
    if (icompdef() != 0) pSpec_N_I = icompdef()->countSpecs();
    if (ocompdef() != 0) pSpec_N_O = ocompdef()->countSpecs();
    
    // 3 -- DEAL WITH PATCH SREAC'S
    if (pSReac_N != 0)
    {
        // Set up local indices.
        pSReac_L2G = new gidxT[countSReacs()];
        uint ngsreacs = statedef()->countSReacs();
        for (uint i = 0; i < ngsreacs; ++i)
        {
            lidxT lidx = sreacG2L(i);
            if (lidx == LIDX_UNDEFINED) continue;
            pSReac_L2G[lidx] = i;
        }
        
        // Create _DEP, _LHS and _UPD vectors.
        uint arrsize_i = 0;
        uint arrsize_s = countSpecs() * countSReacs();
        uint arrsize_o = 0;
        pSReac_DEP_S_Spec = new depT[arrsize_s];
        pSReac_LHS_S_Spec = new uint[arrsize_s];
        pSReac_UPD_S_Spec = new int[arrsize_s];
        // DEBUG: 08-Feb-2008
        std::fill_n(pSReac_DEP_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_LHS_S_Spec, arrsize_s, 0);
        std::fill_n(pSReac_UPD_S_Spec, arrsize_s, 0);
        if (icompdef() != 0) // Only create if inner comp exists.
        {
            arrsize_i = countSpecs_I() * countSReacs();
            pSReac_DEP_I_Spec = new depT[arrsize_i];
            pSReac_LHS_I_Spec = new uint[arrsize_i];
            pSReac_UPD_I_Spec = new int[arrsize_i];
            // DEBUG: 08-Feb-2008
            std::fill_n(pSReac_DEP_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_LHS_I_Spec, arrsize_i, 0);
            std::fill_n(pSReac_UPD_I_Spec, arrsize_i, 0);
        }
        if (ocompdef() != 0) // Only create if outer comp exists.
        {
            arrsize_o = countSpecs_O() * countSReacs();
            pSReac_DEP_O_Spec = new depT[arrsize_o];
            pSReac_LHS_O_Spec = new uint[arrsize_o];
            pSReac_UPD_O_Spec = new int[arrsize_o];
            // DEBUG: 08-Feb-2008
            std::fill_n(pSReac_DEP_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_LHS_O_Spec, arrsize_o, 0);
            std::fill_n(pSReac_UPD_O_Spec, arrsize_o, 0);
        }
        
        // Fill the vectors with all kinds of useful information.
        for (uint ri = 0; ri < countSReacs(); ++ri)
        {
            SReacDef * sreacdef = sreac(ri);
            
            // Handle surface stuff.
            for (uint si = 0; si < ngspecs; ++si)
            {
                if (sreacdef->req_S(si) == false) continue;
                
                // TODO: turn into error check?
                lidxT sil = specG2L(si);
                assert(sil != LIDX_UNDEFINED);
                
                uint aridx = _IDX_SReac_S_Spec(ri, sil);
                pSReac_DEP_S_Spec[aridx] = sreacdef->dep_S(si);
                pSReac_LHS_S_Spec[aridx] = sreacdef->lhs_S(si);
                pSReac_UPD_S_Spec[aridx] = sreacdef->upd_S(si);
            }
            
            // Handle the inside comp stuff.
            if (sreacdef->reqInside() == true)
            {
                // TODO: turn into real error check?
                assert(icompdef() != 0);
                
                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (sreacdef->req_I(si) == false) continue;
                    
                    // TODO: turn into error check?
                    lidxT sil = specG2L_I(si);
                    assert(sil != LIDX_UNDEFINED);
                    
                    uint aridx = _IDX_SReac_I_Spec(ri, sil);
                    pSReac_DEP_I_Spec[aridx] = sreacdef->dep_I(si);
                    pSReac_LHS_I_Spec[aridx] = sreacdef->lhs_I(si);
                    pSReac_UPD_I_Spec[aridx] = sreacdef->upd_I(si);
                }
            }
            
            // Handle the outside comp stuff.
            if (sreacdef->reqOutside() == true)
            {
                // TODO: turn into real error check?
                assert(ocompdef() != 0);
                
                for (uint si = 0; si < ngspecs; ++si)
                {
                    if (sreacdef->req_O(si) == false) continue;
                    
                    // TODO: turn into error check?
                    lidxT sil = specG2L_O(si);
                    assert(sil != LIDX_UNDEFINED);
                    
                    uint aridx = _IDX_SReac_O_Spec(ri, sil);
                    pSReac_DEP_O_Spec[aridx] = sreacdef->dep_O(si);
                    pSReac_LHS_O_Spec[aridx] = sreacdef->lhs_O(si);
                    pSReac_UPD_O_Spec[aridx] = sreacdef->upd_O(si);
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void PatchDef::setupDependencies(void)
{
    // Currently doesn't do anything.
}

////////////////////////////////////////////////////////////////////////////////

SpecDef * PatchDef::spec(lidxT idx) const
{
    assert(pLocalIndicesSetupDone == true);
    return statedef()->spec(specL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

SReacDef * PatchDef::sreac(lidxT idx) const
{
    assert(pLocalIndicesSetupDone == true);
    return statedef()->sreac(sreacL2G(idx));
}

////////////////////////////////////////////////////////////////////////////////

depT PatchDef::sreac_dep_I(lidxT sreac, lidxT spec) const
{
    return pSReac_DEP_I_Spec[spec + (sreac * countSpecs_I())];
}

////////////////////////////////////////////////////////////////////////////////

depT PatchDef::sreac_dep_S(lidxT sreac, lidxT spec) const
{
    return pSReac_DEP_S_Spec[spec + (sreac * countSpecs())];
}

////////////////////////////////////////////////////////////////////////////////

depT PatchDef::sreac_dep_O(lidxT sreac, lidxT spec) const
{
    return pSReac_DEP_O_Spec[spec + (sreac * countSpecs_O())];
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_I_bgn(lidxT sreac) const
{
    return pSReac_LHS_I_Spec + (sreac * countSpecs_I());
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_I_end(lidxT sreac) const
{
    return pSReac_LHS_I_Spec + ((sreac + 1) * countSpecs_I());
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_S_bgn(lidxT sreac) const
{
    return pSReac_LHS_S_Spec + (sreac * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_S_end(lidxT sreac) const
{
    return pSReac_LHS_S_Spec + ((sreac + 1) * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_O_bgn(lidxT sreac) const
{
    return pSReac_LHS_O_Spec + (sreac * countSpecs_O());
}

////////////////////////////////////////////////////////////////////////////////

uint * PatchDef::sreac_lhs_O_end(lidxT sreac) const
{
    return pSReac_LHS_O_Spec + ((sreac + 1) * countSpecs_O());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_I_bgn(lidxT sreac) const
{
    return pSReac_UPD_I_Spec + (sreac * countSpecs_I());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_I_end(lidxT sreac) const
{
    return pSReac_UPD_I_Spec + ((sreac + 1) * countSpecs_I());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_S_bgn(lidxT sreac) const
{
    return pSReac_UPD_S_Spec + (sreac * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_S_end(lidxT sreac) const
{
    return pSReac_UPD_S_Spec + ((sreac + 1) * countSpecs());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_O_bgn(lidxT sreac) const
{
    return pSReac_UPD_O_Spec + (sreac * countSpecs_O());
}

////////////////////////////////////////////////////////////////////////////////

int * PatchDef::sreac_upd_O_end(lidxT sreac) const
{
    return pSReac_UPD_O_Spec + ((sreac + 1) * countSpecs_O());
}

////////////////////////////////////////////////////////////////////////////////

// END
