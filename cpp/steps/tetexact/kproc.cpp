////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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
//
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <vector>
#include <cassert>

// STEPS headers.
#include <steps/common.h>
#include <steps/tetexact/kproc.hpp>

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::tetexact, stex);

////////////////////////////////////////////////////////////////////////////////

stex::KProc::KProc(void)
: rExtent(0)
, pFlags(0)
, pSchedIDX(0)
{
}

////////////////////////////////////////////////////////////////////////////////

stex::KProc::~KProc(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void stex::KProc::setActive(bool active)
{
    if (active == true) pFlags &= ~INACTIVATED;
    else pFlags |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

uint stex::KProc::getExtent(void) const
{
	return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void stex::KProc::resetExtent(void)
{
	rExtent = 0;
}
////////////////////////////////////////////////////////////////////////////////

void stex::KProc::resetCcst(void) const
{
	// This should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

double stex::KProc::c(void) const
{
    // Should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

double stex::KProc::h(void) const
{
    // Should never get called on base object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////
/*
steps::solver::Reacdef * swmd::KProc::defr(void) const
{
	// Should only be called on derived object
	assert (false);
}

////////////////////////////////////////////////////////////////////////////////

steps::solver::SReacdef * swmd::KProc::defsr(void) const
{
	// Should olny be called on derived object
	assert (false);
}
*/
////////////////////////////////////////////////////////////////////////////////

// END
