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

#ifndef STEPS_SIM_SWIGINF_FUNC_CORE_HPP
#define STEPS_SIM_SWIGINF_FUNC_CORE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class State;

////////////////////////////////////////////////////////////////////////////////
// SOLVER INFORMATION
////////////////////////////////////////////////////////////////////////////////

extern char *   siGetSolverName(void);
extern char *   siGetSolverDesc(void);
extern char *   siGetSolverAuthors(void);
extern char *   siGetSolverEmail(void);

////////////////////////////////////////////////////////////////////////////////
// CREATION & DESTRUCTION
////////////////////////////////////////////////////////////////////////////////

extern State *  siNewState(void);
extern void     siDelState(State * s);

extern void     siBeginStateDef(State * s);
extern void     siEndStateDef(State * s);

extern void     siBeginVarDef(State * s);
extern void     siEndVarDef(State * s);
extern uint     siNewSpec(State * s, char * name);

extern void     siBeginReacDef(State * s);
extern void     siEndReacDef(State * s);
extern uint     siNewReac(State * s, char * name, double kcst);
extern void     siAddReacLHS(State * s, uint ridx, uint sidx);
extern void     siAddReacRHS(State * s, uint ridx, uint sidx);

extern void     siBeginDiffDef(State * s);
extern void     siEndDiffDef(State * s);
extern uint     siNewDiff(State * s, char * name, uint sidx, double dcst);

extern void     siBeginSReacDef(State * s);
extern void     siEndSReacDef(State * s);
extern uint     siNewSReac(State * s, char * name, double kcst, bool inside);
extern void     siAddSReacLHS_I(State * s, uint ridx, uint sidx);
extern void     siAddSReacLHS_S(State * s, uint ridx, uint sidx);
extern void     siAddSReacLHS_O(State * s, uint ridx, uint sidx);
extern void     siAddSReacRHS_I(State * s, uint ridx, uint sidx);
extern void     siAddSReacRHS_S(State * s, uint ridx, uint sidx);
extern void     siAddSReacRHS_O(State * s, uint ridx, uint sidx);

extern void     siBeginCompDef(State * s);
extern void     siEndCompDef(State * s);
extern uint     siNewComp(State * s, char * name, double vol);
extern void     siAddCompSpec(State * s, uint cidx, uint sidx);
extern void     siAddCompReac(State * s, uint cidx, uint ridx);
extern void     siAddCompDiff(State * s, uint cidx, uint didx);

extern void     siBeginPatchDef(State * s);
extern void     siEndPatchDef(State * s);
extern uint     siNewPatch(State * s, char * name, double area, 
                           uint cidx_in, uint cidx_out);
extern void     siAddPatchSpec(State * s, uint pidx, uint sidx);
extern void     siAddPatchSReac(State * s, uint pidx, uint ridx);

extern void     siSetRNG(State * s, steps::rng::RNG * rng);

////////////////////////////////////////////////////////////////////////////////
// SIMULATION CONTROLS
////////////////////////////////////////////////////////////////////////////////

extern void     siReset(State * s);
extern void     siRun(State * s, double endtime);
extern void		siSetDT(State * s, double dt);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      GENERAL
////////////////////////////////////////////////////////////////////////////////

extern double   siGetTime(State * s);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      COMPARTMENT
////////////////////////////////////////////////////////////////////////////////

extern double   siGetCompVol(State * s, uint cidx);
extern void     siSetCompVol(State * s, uint cidx, double vol);

extern double   siGetCompCount(State * s, uint cidx, uint sidx);				
extern void     siSetCompCount(State * s, uint cidx, uint sidx, double n);		

extern double   siGetCompAmount(State * s, uint cidx, uint sidx);
extern void     siSetCompAmount(State * s, uint cidx, uint sidx, double m);

extern double   siGetCompConc(State * s, uint cidx, uint sidx);
extern void     siSetCompConc(State * s, uint cidx, uint sidx, double c);

extern bool     siGetCompClamped(State * s, uint cidx, uint sidx);
extern void     siSetCompClamped(State * s, uint cidx, uint sidx, bool buf);

extern double   siGetCompReacK(State * s, uint cidx, uint ridx);
extern void     siSetCompReacK(State * s, uint cidx, uint ridx, double kf);

extern bool     siGetCompReacActive(State * s, uint cidx, uint ridx);
extern void     siSetCompReacActive(State * s, uint cidx, uint ridx, bool act);

extern double 	siGetCompDiffD(State * s, uint cidx, uint didx);
extern void 	siSetCompDiffD(State * s, uint cidx, uint didx);

extern bool		siGetCompDiffActive(State * s, uint cidx, uint didx);
extern void 	siSetCompDiffActive(State * s, uint cidx, uint didx, bool act);

////////////////////////////////////////////////////////////////////////////////
// SOLVER STATE ACCESS: 
//      PATCH
////////////////////////////////////////////////////////////////////////////////

extern double   siGetPatchArea(State * s, uint pidx);
extern void     siSetPatchArea(State * s, uint pidx, double area);

extern double   siGetPatchCount(State * s, uint pidx, uint sidx);
extern void     siSetPatchCount(State * s, uint pidx, uint sidx, double n);

extern double   siGetPatchAmount(State * s, uint pidx, uint sidx);
extern void     siSetPatchAmount(State * s, uint pidx, uint sidx, double m);

extern bool     siGetPatchClamped(State * s, uint pidx, uint sidx);
extern void     siSetPatchClamped(State * s, uint pidx, uint sidx, bool buf);

extern double   siGetPatchSReacK(State * s, uint pidx, uint ridx);
extern void     siSetPatchSReacK(State * s, uint pidx, uint ridx, double kf);

extern bool     siGetPatchSReacActive(State * s, uint pidx, uint ridx);
extern void     siSetPatchSReacActive(State * s, uint pidx, uint ridx, bool a);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SWIGINF_FUNC_CORE_HPP

// END
