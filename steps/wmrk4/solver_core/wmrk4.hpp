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
// 
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WMRK4_SOLVER_CORE_WMRK4_HPP
#define STEPS_WMRK4_SOLVER_CORE_WMRK4_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/constants.hpp>
#include <steps/wmrk4/solver_core/state.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class State;

// Auxiliary declarations.
typedef std::vector<double>             dVec;
typedef dVec::iterator					dVecI;

typedef std::vector<uint>				uiVec;
typedef uiVec::iterator					uiVecI;

////////////////////////////////////////////////////////////////////////////////

class Wmrk4
	{
		
	public: 
		
		////////////////////////////////////////////////////////////////////////
		
		Wmrk4(State * s);														
		~Wmrk4(void);
		
		////////////////////////////////////////////////////////////////////////
		
		/// initialises vectors and builds the reaction matrix. 
		///
		void setup(void);
		
		/// this function refills the values and flags vectors 
		/// called if a flag or a count changed
		///
		void refill(void);
		
		/// returns properly scaled reaction constant
		///
		double comp_ccst(double kcst, double vol, uint order);
		
		/// the Runge-Kutta algorithm
		///
		void rk4(void);							
		
		/// the simple stepper
		///
		void rksteps(double t1, double t2);		
		
		/// the derivatives calculator
		///
		void setderivs(dVec& vals, dVec& dydx);
		
		/// update local values vector, 
		/// then update state with computed counts
		///
		void update(void);
		
		////////////////////////////////////////////////////////////////////////
		
		State * state(void) const
		{ return pState; }
		////////////////////////////////////////////////////////////////////////
		
		
		
	private:
		
		////////////////////////////////////////////////////////////////////////	
		
		State *						pState;
		
		/// the reaction matrix 
		uint **						pReacMtx;					
		
		/// update matrix of rhs - lhs
		int **						pUpdMtx;				
		
		/// number of species total: all species in all comps and patches
		uint						pSpecs_tot;					
		
		/// number of reactions total: all reactions and surface reactions
		/// in each comp and patch
		uint						pReacs_tot;	
		
		/// Properly scaled reaction constant vector
		dVec						pCcst;	
		
		/// vector holding current molecular counts (as doubles)
		dVec						pVals;						
		
		/// vector holding flags on species
		uiVec						pSFlags;
		
		/// vector holding flags on reactions
		uiVec						pRFlags;
		
		/// vector holding new, calculated counts
		dVec						pNewVals;					
		
		/// vector of present derivatives
		dVec						pDyDx;						
		
		/// matrix of derivatives from each reaction for each species
		double **					pDyDxlhs;					
		
		/// the time step 
		double						pdt;
		
		/// objects to contain temporary values important in algorithm
		dVec						yt;
		dVec						dyt;
		dVec						dym;
		////////////////////////////////////////////////////////////////////////	
		
	};	

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_WMRK4_SOLVER_CORE_WMRK4_HPP

// END


