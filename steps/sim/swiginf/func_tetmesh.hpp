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

#ifndef STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP
#define STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP 1

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
// CREATION & DESTRUCTION
////////////////////////////////////////////////////////////////////////////////

void siBeginTetmeshDef(State * s);
void siEndTetmeshDef(State *s);

void siBeginTetDef(State * s);
void siEndTetDef(State * s);
uint siNewTet(State * s, uint cidx, double vol, 
	double a1, double a2, double a3, double a4,
	double d1, double d2, double d3, double d4);

void siBeginConnectDef(State * s);
void siEndConnectDef(State * s);
void siConnectTetTet(State * s, uint side, uint tidx1, uint tidx2);
void siConnectTetTriInside(State * s, uint side, uint tetidx, uint triidx);
void siConnectTetTriOutside(State * s, uint side, uint tetidx, uint triidx);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      TETRAHEDRAL VOLUME ELEMENTS
////////////////////////////////////////////////////////////////////////////////

extern double   siGetTetVol(State * s, uint tidx);
extern void     siSetTetVol(State * s, uint tidx, double vol);

extern uint     siGetTetCount(State * s, uint tidx, uint sidx);
extern void     siSetTetCount(State * s, uint tidx, uint sidx, uint n);

extern double   siGetTetMass(State * s, uint tidx, uint sidx);
extern void     siSetTetMass(State * s, uint tidx, uint sidx, double m);

extern double   siGetTetConc(State * s, uint tidx, uint sidx);
extern void     siSetTetConc(State * s, uint tidx, uint sidx, double c);

extern bool     siGetTetClamped(State * s, uint tidx, uint sidx);
extern void     siSetTetClamped(State * s, uint tidx, uint sidx, bool buf);

extern double   siGetTetReacK(State * s, uint tidx, uint ridx);
extern void     siSetTetReacK(State * s, uint tidx, uint ridx, double kf);

extern bool     siGetTetReacActive(State * s, uint tidx, uint ridx);
extern void     siSetTetReacActive(State * s, uint tidx, uint ridx, bool act);

extern double 	siGetTetDiffD(State * s, uint tidx, uint didx);
extern void 	siSetTetDiffD(State * s, uint tidx, uint didx);

extern bool		siGetTetDiffActive(State * s, uint tidx, uint didx);
extern void 	siGetTetDiffActive(State * s, uint tidx, uint didx, bool act);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP

// END
