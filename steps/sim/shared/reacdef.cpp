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
#include <steps/math/accumulate.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/statedef.hpp>

NAMESPACE_ALIAS(steps::math, smath);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

ReacDef::ReacDef(StateDef * sdef, gidxT idx, std::string const & name)
: pStateDef(sdef)
, pGIDX(idx)
, pName(name)
, pOrder(0)
, pKcst(0.0)
, pFinalSetupDone(false)
, pSpec_DEP(0)
, pSpec_LHS(0)
, pSpec_RHS(0)
, pSpec_UPD(0)
, pSpec_UPD_Coll()
{
    uint nspecs = statedef()->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_DEP = new depT[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);
    pSpec_LHS = new uint[nspecs];
    std::fill_n(pSpec_LHS, nspecs, 0);
    pSpec_RHS = new uint[nspecs];
    std::fill_n(pSpec_RHS, nspecs, 0);
    pSpec_UPD = new int[nspecs];
    std::fill_n(pSpec_UPD, nspecs, 0);
}

////////////////////////////////////////////////////////////////////////////////

ReacDef::~ReacDef(void)
{
    delete[] pSpec_DEP;
    delete[] pSpec_LHS;
    delete[] pSpec_RHS;
    delete[] pSpec_UPD;
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::incLHS(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_LHS[idx] += 1;
    computeOrder();
    return pSpec_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setLHS(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_LHS[idx] = n;
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::incRHS(gidxT idx)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_RHS[idx] += 1;
    return pSpec_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setRHS(gidxT idx, uint n)
{
    assert(pFinalSetupDone == false);
    assert(idx < statedef()->countSpecs());
    pSpec_RHS[idx] = n;
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setupFinal(void)
{
    assert(pFinalSetupDone == false);
    
    computeOrder();
    uint nspecs = statedef()->countSpecs();
    for (uint i = 0; i < nspecs; ++i)
    {
        int lhs = static_cast<int>(pSpec_LHS[i]); 
        int aux = pSpec_UPD[i] = static_cast<int>(pSpec_RHS[i]) - lhs;
        if (lhs != 0) pSpec_DEP[i] |= DEP_STOICH;
        if (aux != 0) pSpec_UPD_Coll.push_back(i);
    }
    
    pFinalSetupDone = true;
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::lhs(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_LHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

depT ReacDef::dep(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_DEP[idx];
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::rhs(gidxT idx) const
{
    assert(idx < statedef()->countSpecs());
    return pSpec_RHS[idx];
}

////////////////////////////////////////////////////////////////////////////////

int ReacDef::upd(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_UPD[idx];
}

////////////////////////////////////////////////////////////////////////////////

bool ReacDef::req(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (pSpec_DEP[idx] != DEP_NONE) return true;
    if (pSpec_RHS[idx] != 0) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::computeOrder(void)
{
    // Compute the order of the reaction.
    uint nspecs = statedef()->countSpecs();
    pOrder = smath::accumulate(pSpec_LHS, pSpec_LHS + nspecs, 0);
}

////////////////////////////////////////////////////////////////////////////////

// END
