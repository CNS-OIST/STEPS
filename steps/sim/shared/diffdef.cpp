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
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/shared/types.hpp>

USING(std, string);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

DiffDef::DiffDef(StateDef * sdef, gidxT idx, string const & name)
: pStateDef(sdef)
, pFinalSetupDone(false)
, pGIDX(idx)
, pName(name)
, pDcst(0.0)
, pSpec_DEP(0)
, pSpec_LIG(GIDX_UNDEFINED)
{
    uint nspecs = statedef()->countSpecs();
    if (nspecs == 0) return; // Would be weird, but okay.
    pSpec_DEP = new depT[nspecs];
    std::fill_n(pSpec_DEP, nspecs, DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

DiffDef::~DiffDef(void)
{
    delete[] pSpec_DEP;
}

////////////////////////////////////////////////////////////////////////////////

void DiffDef::setLig(gidxT idx)
{
    assert(pFinalSetupDone == false);
    pSpec_LIG = idx;
}

////////////////////////////////////////////////////////////////////////////////

void DiffDef::setupFinal(void)
{
    assert(pFinalSetupDone == false);
    pSpec_DEP[lig()] = DEP_STOICH;
    pFinalSetupDone = true;
}

////////////////////////////////////////////////////////////////////////////////

void DiffDef::setDcst(double const & d)
{
    assert(pFinalSetupDone == false);
    pDcst = d;
}

////////////////////////////////////////////////////////////////////////////////

depT DiffDef::dep(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    return pSpec_DEP[idx];
}

////////////////////////////////////////////////////////////////////////////////

bool DiffDef::req(gidxT idx) const
{
    assert(pFinalSetupDone == true);
    assert(idx < statedef()->countSpecs());
    if (pSpec_DEP[idx] != DEP_NONE) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////
/*
bool DiffDef::dependsOnSpec(uint gidx) const
{
    assert(gidx < statedef()->countSpecs());
    if (gidx != pLig) return false;
    return true;
}
*/
////////////////////////////////////////////////////////////////////////////////
/*
bool DiffDef::affectsSpec(uint gidx) const
{
    assert(gidx < statedef()->countSpecs());
    if (gidx != pLig) return false;
    return true;
}
*/
////////////////////////////////////////////////////////////////////////////////

// END
