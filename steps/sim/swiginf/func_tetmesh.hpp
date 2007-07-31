////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP
#define STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP 1

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
// TODO: improve interface here! E.g. why do we need to specify areas and
// distances already here, when we have to specificy connectivity later?
uint siNewTet(State * s, uint cidx, double vol, 
	double a1, double a2, double a3, double a4,
	double d1, double d2, double d3, double d4);

// TODO: rename the following functions
void siBeginConnectDef(State * s);
void siEndConnectDef(State * s);
// TODO: make it possible to declare tet-tet neighbourship at once
// for a given tetrahedron.
void siConnectTetTet(State * s, uint side, uint tidx1, uint tidx2);
void siConnectTetTriInside(State * s, uint side, uint tetidx, uint triidx);
void siConnectTetTriOutside(State * s, uint side, uint tetidx, uint triidx);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SWIGINF_FUNC_TETMESH_HPP

// END
