////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_SWIGINF_FUNC_CORE_HPP
#define STEPS_SIM_SWIGINF_FUNC_CORE_HPP 1

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

extern void     siBeginCompDef(State * s);
extern void     siEndCompDef(State * s);
extern uint     siNewComp(State * s, char * name, double vol);
extern void     siAddCompSpec(State * s, uint cidx, uint sidx);
extern void     siAddCompReac(State * s, uint cidx, uint ridx);
extern void     siAddCompDiff(State * s, uint cidx, uint didx);

extern void     siSetRNG(State * s, steps::rng::RNG * rng);


////////////////////////////////////////////////////////////////////////////////
// SIMULATION CONTROLS
////////////////////////////////////////////////////////////////////////////////

extern void     siReset(State * s);
extern void     siRun(State * s, double endtime);

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

extern uint     siGetCompCount(State * s, uint cidx, uint sidx);
extern void     siSetCompCount(State * s, uint cidx, uint sidx, uint n);

extern double   siGetCompMass(State * s, uint cidx, uint sidx);
extern void     siSetCompMass(State * s, uint cidx, uint sidx, double m);

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
extern void 	siGetCompDiffActive(State * s, uint cidx, uint didx, bool act);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_SWIGINF_FUNC_CORE_HPP

// END
