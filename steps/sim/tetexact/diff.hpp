////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
// 
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_SIM_TETEXACT_DIFF_HPP
#define STEPS_SIM_TETEXACT_DIFF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/tetexact/kproc.hpp>
#include <steps/sim/tetexact/sched.hpp>

// Forward declarations.
class DiffDef;
class State;
class Tet;

////////////////////////////////////////////////////////////////////////////////

class Diff
: public KProc
{

public:
	
	////////////////////////////////////////////////////////////////////////
	
	/// Constructor.
	///
	Diff(DiffDef * ddef, Tet * tet);
	
	/// Destructor.
	///
	virtual ~Diff(void);
	
	///
	virtual void setupDeps(void);
	
	virtual bool depSpecTet(uint gidx, Tet * tet);
	
	virtual void reset(void);
	
	////////////////////////////////////////////////////////////////////////
	
	///
	virtual double rate(void) const;
	
	///
	virtual SchedIDXVec const & apply(State * s);

	////////////////////////////////////////////////////////////////////////
	
	inline DiffDef * def(void) const
	{ return pDiffDef; }
	
	////////////////////////////////////////////////////////////////////////
	
private:

	////////////////////////////////////////////////////////////////////////
	
	DiffDef * 					pDiffDef;
	
	Tet * 						pTet;
	
	/// 
	SchedIDXVec              	pUpdVec[4];
    
    ////////////////////////////////////////////////////////////////////////
    
	/// Properly scaled diffusivity constant.    
    ///
	double                      pScaledDcst;
    
    double                      pCDFSelector[3];	
};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SIM_TETEXACT_DIFF_HPP

// END
