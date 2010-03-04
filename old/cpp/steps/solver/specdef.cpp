////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2010ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/solver/types.hpp>
#include <steps/error.hpp>
#include <steps/solver/statedef.hpp>
#include <steps/solver/specdef.hpp>
#include <steps/model/spec.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::solver, ssolver);

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::Specdef(Statedef * sd, uint idx, steps::model::Spec * s)
: pStatedef(sd)
, pIdx(idx)
, pName()
, pSetupdone(false)
{
	assert(pStatedef != 0);
	assert(s != 0);
	pName = s->getID();

}

////////////////////////////////////////////////////////////////////////////////

ssolver::Specdef::~Specdef(void)
{

}

////////////////////////////////////////////////////////////////////////////////

std::string const ssolver::Specdef::name(void) const
{
	return pName;
}

////////////////////////////////////////////////////////////////////////////////

void ssolver::Specdef::setup(void)
{

}

////////////////////////////////////////////////////////////////////////////////


// END

