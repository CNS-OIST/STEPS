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
#include <cassert>
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/accumulate.hpp>
#include <steps/sim/shared/sreacdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>

NAMESPACE_ALIAS(steps::math, smath);
USING(std, string);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

SReacDef::SReacDef(StateDef * sdef, gidxT idx, string const & name, orientT o)
: pStateDef(sdef)
, pGIDX(idx)
, pName(name)
, pOrient(o)
, pOrder(0)
, pKcst(0.0)
, pFinalSetupDone(false)
, pSpec_I_DEP(0)
, pSpec_S_DEP(0)
, pSpec_O_DEP(0)
, pSpec_I_LHS(0)
, pSpec_S_LHS(0)
, pSpec_O_LHS(0)
, pSpec_I_RHS(0)
, pSpec_S_RHS(0)
, pSpec_O_RHS(0)
, pSpec_I_UPD(0)
, pSpec_S_UPD(0)
, pSpec_O_UPD(0)
, pSpec_I_UPD_Coll()
, pSpec_S_UPD_Coll()
, pSpec_O_UPD_Coll()
{
    uint nspecs = statedef()->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_S_DEP = new depT[nspecs];
    std::fill_n(pSpec_S_DEP, nspecs, DEP_NONE);
    pSpec_S_LHS = new uint[nspecs];
    std::fill_n(pSpec_S_LHS, nspecs, 0);
    if (pOrient == SReacDef::INSIDE)
    {
        pSpec_I_DEP = new depT[nspecs];
        std::fill_n(pSpec_I_DEP, nspecs, DEP_NONE);
        pSpec_I_LHS = new uint[nspecs];
        std::fill_n(pSpec_I_LHS, nspecs, 0);
    }
    else
    {
        pSpec_O_DEP = new depT[nspecs];
        std::fill_n(pSpec_O_DEP, nspecs, DEP_NONE);
        pSpec_O_LHS = new uint[nspecs];
        std::fill_n(pSpec_O_LHS, nspecs, 0);        
    }
    pSpec_I_RHS = new uint[nspecs];
    std::fill_n(pSpec_I_RHS, nspecs, 0);
    pSpec_S_RHS = new uint[nspecs];
    std::fill_n(pSpec_S_RHS, nspecs, 0);
    pSpec_O_RHS = new uint[nspecs];
    std::fill_n(pSpec_O_RHS, nspecs, 0);
    pSpec_I_UPD = new int[nspecs];
    std::fill_n(pSpec_I_UPD, nspecs, 0);
    pSpec_S_UPD = new int[nspecs];
    std::fill_n(pSpec_S_UPD, nspecs, 0);
    pSpec_O_UPD = new int[nspecs];
    std::fill_n(pSpec_O_UPD, nspecs, 0);    
}

////////////////////////////////////////////////////////////////////////////////

SReacDef::~SReacDef(void)
{
    delete[] pSpec_I_DEP;
    delete[] pSpec_S_DEP;
    delete[] pSpec_O_DEP;
    delete[] pSpec_I_LHS;
    delete[] pSpec_S_LHS;
    delete[] pSpec_O_LHS;
    delete[] pSpec_I_RHS;
    delete[] pSpec_S_RHS;
    delete[] pSpec_O_RHS;
    delete[] pSpec_I_UPD;
    delete[] pSpec_S_UPD;
    delete[] pSpec_O_UPD;
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incLHS_I(gidxT idx)
{
    assert(pFinalSetupDone == false);
    if (outside()) return 0;
    assert(idx < statedef()->countSpecs());
    pSpec_I_LHS[idx] += 1;
    computeOrder();
    return pSpec_I_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incLHS_S(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_S_LHS[idx] += 1;
    computeOrder();
    return pSpec_S_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incLHS_O(gidxT idx)
{
    assert(pFinalSetupDone == false);
    if (inside()) return 0;
    assert(idx < statedef()->countSpecs());
    pSpec_O_LHS[idx] += 1;
    computeOrder();
    return pSpec_O_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setLHS_I(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    if (outside()) return;
    assert(idx < statedef()->countSpecs());
    pSpec_I_LHS[idx] = n;
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setLHS_S(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_S_LHS[idx] = n;
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setLHS_O(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    if (inside()) return;
    assert(idx < statedef()->countSpecs());
    pSpec_O_LHS[idx] = n;
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incRHS_I(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_I_RHS[idx] += 1;
    return pSpec_I_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incRHS_S(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_S_RHS[idx] += 1;
    return pSpec_S_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::incRHS_O(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_O_RHS[idx] += 1;
    return pSpec_O_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setRHS_I(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_I_RHS[idx] = n;
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setRHS_S(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_S_RHS[idx] = n;
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setRHS_O(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_O_RHS[idx] = n;
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setupFinal(void)
{
    assert(pFinalSetupDone == false);
    
    computeOrder();
    uint nspecs = statedef()->countSpecs();
    
    // Deal with surface.
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = static_cast<int>(pSpec_S_LHS[i]); 
        int aux = pSpec_S_UPD[i] = static_cast<int>(pSpec_S_RHS[i]) - lhs;
        if (lhs != 0) pSpec_S_DEP[i] |= DEP_STOICH;
        if (aux != 0) pSpec_S_UPD_Coll.push_back(i);
    }
    
    // Deal with inside.
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = (inside() ? static_cast<int>(pSpec_I_LHS[i]) : 0); 
        int aux = pSpec_I_UPD[i] = static_cast<int>(pSpec_I_RHS[i]) - lhs;
        if (lhs != 0) pSpec_I_DEP[i] |= DEP_STOICH;
        if (aux != 0) pSpec_I_UPD_Coll.push_back(i);
    }
    
    // Deal with outside.
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = (outside() ? static_cast<int>(pSpec_O_LHS[i]) : 0); 
        int aux = pSpec_O_UPD[i] = static_cast<int>(pSpec_O_RHS[i]) - lhs;
        if (lhs != 0) pSpec_O_DEP[i] |= DEP_STOICH;
        if (aux != 0) pSpec_O_UPD_Coll.push_back(i);
    }
    
    pFinalSetupDone = true;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacDef::reqInside(void) const
{
    assert(pFinalSetupDone == true);
    
    // This can be checked by seeing if DEP_I or RHS_I is non-zero 
    // for any species.
    uint nspecs = statedef()->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
        if (req_I(i) == true) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacDef::reqOutside(void) const
{
    assert(pFinalSetupDone == true);
    
    // This can be checked by seeing if DEP_I or RHS_I is non-zero 
    // for any species.
    uint nspecs = statedef()->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
        if (req_O(i) == true) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

double SReacDef::kcst(void) const
{
    return pKcst;
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::setKcst(double const & k)
{
    assert(pFinalSetupDone == false);
    pKcst = k;
}
    
////////////////////////////////////////////////////////////////////////////////

uint SReacDef::lhs_I(gidxT idx) const
{
    if (outside()) return 0;
    assert(idx < statedef()->countSpecs());
    return pSpec_I_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::lhs_S(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_S_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::lhs_O(gidxT idx) const
{
    if (inside()) return 0;
    assert(idx < statedef()->countSpecs());
    return pSpec_O_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

depT SReacDef::dep_I(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (outside()) return DEP_NONE;
    return pSpec_I_DEP[idx]; 
}

////////////////////////////////////////////////////////////////////////////////

depT SReacDef::dep_S(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_S_DEP[idx];
}

////////////////////////////////////////////////////////////////////////////////

depT SReacDef::dep_O(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (inside()) return DEP_NONE;
    return pSpec_O_DEP[idx];
}
    
////////////////////////////////////////////////////////////////////////////////

uint SReacDef::rhs_I(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_I_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::rhs_S(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_S_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint SReacDef::rhs_O(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_O_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

int SReacDef::upd_I(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_I_UPD[idx];
}

////////////////////////////////////////////////////////////////////////////////

int SReacDef::upd_S(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_S_UPD[idx];
}

////////////////////////////////////////////////////////////////////////////////

int SReacDef::upd_O(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_O_UPD[idx];
}

////////////////////////////////////////////////////////////////////////////////

bool SReacDef::req_I(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (inside())
        if (pSpec_I_DEP[idx] != DEP_NONE) return true;
    if (pSpec_I_RHS[idx] != 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacDef::req_S(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (pSpec_S_DEP[idx] != DEP_NONE) return true;
    if (pSpec_S_RHS[idx] != 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool SReacDef::req_O(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (outside())
        if (pSpec_O_DEP[idx] != DEP_NONE) return true;
    if (pSpec_O_RHS[idx] != 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void SReacDef::computeOrder(void)
{
    uint nspecs = statedef()->countSpecs();
    pOrder = smath::accumulate(pSpec_S_LHS, pSpec_S_LHS + nspecs, 0);
    if (inside())
    {
        pOrder = smath::accumulate(pSpec_I_LHS, pSpec_I_LHS + nspecs, pOrder);
    }
    else if (outside())
    {
        pOrder = smath::accumulate(pSpec_O_LHS, pSpec_O_LHS + nspecs, pOrder);
    }
    else
    {
        assert(0); // shouldn't happen.
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
