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
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/accumulate.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/statedef.hpp>

NAMESPACE_ALIAS(steps::math, smath);

////////////////////////////////////////////////////////////////////////////////

ReacDef::ReacDef(StateDef * sdef, uint gidx, std::string const & name)
: pStateDef(sdef)
, pGIDX(gidx)
, pName(name)
, pOrder(0)
, pKcst(0.0)
, pLHS(sdef->countSpecs())
, pRHS(sdef->countSpecs())
{
    std::fill(pLHS.begin(), pLHS.end(), 0);
    std::fill(pRHS.begin(), pRHS.end(), 0);
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setupFinal(void)
{
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::lhs(uint gidx) const
{
    return pLHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::incLHS(uint gidx)
{
    pLHS[gidx] += 1;
    computeOrder();
    return pLHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setLHS(uint gidx, uint n)
{
    pLHS[gidx] = n;
    computeOrder();
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::rhs(uint gidx) const
{
    return pRHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

uint ReacDef::incRHS(uint gidx)
{
    pRHS[gidx] += 1;
    return pRHS[gidx];
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::setRHS(uint gidx, uint n)
{
    pRHS[gidx] = n;
}

////////////////////////////////////////////////////////////////////////////////

bool ReacDef::dependsOnSpec(uint gidx) const
{
    assert(gidx < statedef()->countSpecs());
    if (lhs(gidx) == 0) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

bool ReacDef::affectsSpec(uint gidx) const
{
    assert(gidx < statedef()->countSpecs());
    if ((rhs(gidx) - lhs(gidx)) == 0) return false;
    return true;
}

////////////////////////////////////////////////////////////////////////////////

void ReacDef::computeOrder(void)
{
    // Compute the order of the reaction.
    pOrder = smath::accumulate(pLHS.begin(), pLHS.end(), 0);
}

////////////////////////////////////////////////////////////////////////////////

// END
