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
// $Id:kproc.cpp 64 2007-08-20 06:25:41Z stefan $
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/tetexact/solver_core/kproc.hpp>
#include <steps/tetexact/solver_core/sched.hpp>

////////////////////////////////////////////////////////////////////////////////

KProc::KProc(void)
: pFlags(0)
, pSchedIDX(0)
{
}

////////////////////////////////////////////////////////////////////////////////

KProc::~KProc(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void KProc::setActive(bool active)
{
    if (active == true) pFlags &= ~INACTIVATED;
    else pFlags |= INACTIVATED;
}

////////////////////////////////////////////////////////////////////////////////

// END
