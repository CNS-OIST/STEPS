////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
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


/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

// Standard library & STL headers.
// #include <vector>

// STEPS headers.
#include "../common.h"
#include "kproc.hpp"

////////////////////////////////////////////////////////////////////////////////

NAMESPACE_ALIAS(steps::wmdirect, swmd);

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::KProc(void)
: rExtent(0)
// , pFlags(0)
, pSchedIDX(0)
{

}

////////////////////////////////////////////////////////////////////////////////

swmd::KProc::~KProc(void)
{

}

////////////////////////////////////////////////////////////////////////////////

uint swmd::KProc::getExtent(void) const
{
	return rExtent;
}

////////////////////////////////////////////////////////////////////////////////

void swmd::KProc::resetExtent(void)
{
	rExtent = 0;
}

////////////////////////////////////////////////////////////////////////////////

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

////////////////////////////////////////////////////////////////////////////////

// END

