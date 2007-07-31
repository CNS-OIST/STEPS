////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <algorithm>
#include <cassert>
#include <functional>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/tetexact/sched.hpp>
#include <steps/sim/tetexact/state.hpp>
#include <steps/sim/tetexact/tet.hpp>

////////////////////////////////////////////////////////////////////////////////
	
State::State(void)
: pStateDef(0)
, pRNG(0)
, pSched(0)
, pTime(0.0)
, pTets()
{
	pStateDef = new StateDef();
	pSched = new Sched();
}

////////////////////////////////////////////////////////////////////////////////

State::~State(void)
{
	delete pStateDef;
}

////////////////////////////////////////////////////////////////////////////////

void State::setupState(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void State::setupTetmesh(void)
{
	//std::for_each(pTets.begin(), pTets.end(), 
	//	std::bind2nd(std::mem_fun(Tet::setupKProcs), pSched));
	for (std::vector<Tet*>::const_iterator i = pTets.begin(); 
		i != pTets.end(); ++i)
	{
		(*i)->setupKProcs(pSched);
	}
}

////////////////////////////////////////////////////////////////////////////////

void State::reset(void)
{
	std::for_each(pTets.begin(), pTets.end(), std::mem_fun(&Tet::reset));
	pSched->reset();
}

////////////////////////////////////////////////////////////////////////////////

uint State::addTet
(
	CompDef * cdef, double vol, 
	double a1, double a2, double a3, double a4,
	double d1, double d2, double d3, double d4
)
{
	Tet * t = new Tet(cdef, vol, a1, a2, a3, a4, d1, d2, d3, d4);
	uint tidx = pTets.size();
	pTets.push_back(t);
	return tidx;
}

////////////////////////////////////////////////////////////////////////////////

// END
