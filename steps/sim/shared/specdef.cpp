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
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/diffdef.hpp>
#include <steps/sim/shared/reacdef.hpp>
#include <steps/sim/shared/specdef.hpp>
#include <steps/sim/shared/statedef.hpp>

////////////////////////////////////////////////////////////////////////////////

SpecDef::SpecDef(StateDef * sdef, uint gidx, std::string const & name)
: pStateDef(sdef)
, pGIDX(gidx)
, pName(name)
{
}

////////////////////////////////////////////////////////////////////////////////

SpecDef::~SpecDef(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void SpecDef::setupFinal(void)
{
}

////////////////////////////////////////////////////////////////////////////////

bool SpecDef::dependsOnReac(uint gidx) const
{
    assert(gidx < statedef()->countReacs());
    ReacDef * r = statedef()->reac(gidx);
    return r->affectsSpec(this->gidx());
}

////////////////////////////////////////////////////////////////////////////////

bool SpecDef::affectsReac(uint gidx) const
{
    assert(gidx < statedef()->countReacs());
    ReacDef * r = statedef()->reac(gidx);
    return r->dependsOnSpec(this->gidx());
}

////////////////////////////////////////////////////////////////////////////////

bool SpecDef::dependsOnDiff(uint gidx) const
{
    assert(gidx < statedef()->countDiffs());
    DiffDef * d = statedef()->diff(gidx);
    return d->affectsSpec(this->gidx());
}

////////////////////////////////////////////////////////////////////////////////

bool SpecDef::affectsDiff(uint gidx) const
{
    assert(gidx < statedef()->countDiffs());
    DiffDef * d = statedef()->diff(gidx);
    return d->dependsOnSpec(this->gidx());
}

////////////////////////////////////////////////////////////////////////////////

// END
