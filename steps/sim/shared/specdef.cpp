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
#include <steps/sim/shared/types.hpp>

USING(std, string);
USING_NAMESPACE(steps::sim);

////////////////////////////////////////////////////////////////////////////////

SpecDef::SpecDef(StateDef * sdef, gidxT idx, std::string const & name)
: pStateDef(sdef)
, pGIDX(idx)
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
/*
bool SpecDef::dependsOnReac(gidxT idx) const
{
    assert(idx < statedef()->countReacs());
    ReacDef * r = statedef()->reac(idx);
    return r->affectsSpec(this->idx());
}
*/
////////////////////////////////////////////////////////////////////////////////
/*
bool SpecDef::affectsReac(gidxT idx) const
{
    assert(idx < statedef()->countReacs());
    ReacDef * r = statedef()->reac(idx);
    return r->dependsOnSpec(this->idx());
}
*/
////////////////////////////////////////////////////////////////////////////////
/*
bool SpecDef::dependsOnDiff(gidxT idx) const
{
    assert(idx < statedef()->countDiffs());
    DiffDef * d = statedef()->diff(idx);
    return d->affectsSpec(this->idx());
}
*/
////////////////////////////////////////////////////////////////////////////////
/*
bool SpecDef::affectsDiff(gidxT idx) const
{
    assert(idx < statedef()->countDiffs());
    DiffDef * d = statedef()->diff(idx);
    return d->dependsOnSpec(this->idx());
}
*/
////////////////////////////////////////////////////////////////////////////////

// END
