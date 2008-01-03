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

#ifndef STEPS_SIM_SWIGINF_FUNC_SSA_HPP
#define STEPS_SIM_SWIGINF_FUNC_SSA_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class State;

////////////////////////////////////////////////////////////////////////////////
// SIMULATION CONTROLS
////////////////////////////////////////////////////////////////////////////////

extern double   siStep(State * s);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      GENERAL
////////////////////////////////////////////////////////////////////////////////

extern uint     siGetNSteps(State * s);

extern double   siGetA0(State * s);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      COMPARTMENT
////////////////////////////////////////////////////////////////////////////////

extern double   siGetCompReacC(State * s, uint cidx, uint ridx);
extern double   siGetCompReacH(State * s, uint cidx, uint ridx);
extern double   siGetCompReacA(State * s, uint cidx, uint ridx);

extern uint     siGetCompReacExtent(State * s, uint cidx, uint ridx);
extern void     siResetCompReacExtent(State * s, uint cidx, uint ridx);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      PATCH
////////////////////////////////////////////////////////////////////////////////

extern double   siGetPatchSReacC(State * s, uint pidx, uint ridx);
extern double   siGetPatchSReacH(State * s, uint pidx, uint ridx);
extern double   siGetPatchCReacA(State * s, uint pidx, uint ridx);

extern uint     siGetPatchSReacExtent(State * s, uint pidx, uint ridx);
extern void     siResetPatchSReacExtent(State * s, uint pidx, uint ridx);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SWIGINF_FUNC_SSA_HPP

// END
