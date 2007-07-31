////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_TETEXACT_STATE_HPP
#define STEPS_SIM_TETEXACT_STATE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/sim/tetexact/sched.hpp>
#include <steps/sim/tetexact/tet.hpp>

////////////////////////////////////////////////////////////////////////////////

class State
{
	
public:
	
	////////////////////////////////////////////////////////////////////////
	
	/// Default constructor.
	///
	State(void);
	
	/// Destructor.
	///
	~State(void);
	
	/// Return the state definition.
	///
	inline StateDef * def(void) const
	{ return pStateDef; }
	
	////////////////////////////////////////////////////////////////////////
	
	void setupState(void);
	
	void setupTetmesh(void);
	
	void reset(void);
	
	////////////////////////////////////////////////////////////////////////
	
	inline void setRNG(steps::rng::RNG * rng)
	{ pRNG = rng; }
	
	inline steps::rng::RNG * rng(void) const
	{ return pRNG; }
	
	////////////////////////////////////////////////////////////////////////
	
	inline double time(void) const
	{ return pTime; }
	
	inline Sched * sched(void) const
	{ return pSched; }
	
	////////////////////////////////////////////////////////////////////////
	
	uint addTet
	(
		CompDef * cdef, double vol, 
		double a1, double a2, double a3, double a4,
		double d1, double d2, double d3, double d4
	);
	
	inline Tet * tet(uint tidx) const
	{ return pTets[tidx]; }
	
private:
	
	StateDef *                  pStateDef;
	
	steps::rng::RNG *           pRNG;
    
	double						pTime;

	Sched * 					pSched;
	
	std::vector<Tet *>			pTets;
	
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_TETEXACT_STATE_HPP

// END
