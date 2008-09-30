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

#ifndef STEPS_WMRK4_SOLVER_CORE_STATE_HPP
#define STEPS_WMRK4_SOLVER_CORE_STATE_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>
#include <steps/sim/shared/compdef.hpp>
#include <steps/sim/shared/patchdef.hpp>
#include <steps/sim/shared/statedef.hpp>
#include <steps/wmrk4/solver_core/comp.hpp>
#include <steps/wmrk4/solver_core/patch.hpp>
#include <steps/wmrk4/solver_core/wmrk4.hpp>										

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Wmrk4;

////////////////////////////////////////////////////////////////////////////////

class State													 
	{															
		
	public:
		
		enum PoolFlags
		{
			CLAMPED_POOLFLAG       = 1
		};
		static const uint PoolFlagDefault = 0;
		
		enum ReacFlags
		{
			INACTIVE_REACFLAG         = 1,
			KCONST_REACFLAG         = 2
		};
		static const uint ReacFlagDefault = INACTIVE_REACFLAG | KCONST_REACFLAG;
		
		////////////////////////////////////////////////////////////////////////
		
		State(void);
		~State(void);
		
		inline steps::sim::StateDef * def(void) const
		{ return pStateDef; }
		
		
		////////////////////////////////////////////////////////////////////////
		
		void setupState(void);
		
		void reset(void);
		
		////////////////////////////////////////////////////////////////////////
		
		void step(void);
		
		void run(double maxt);								
		
		////////////////////////////////////////////////////////////////////////
		
		inline void setRNG(steps::rng::RNG * rng)
		{ pRNG = rng; }
		
		inline steps::rng::RNG * rng(void) const
		{ return pRNG; }
		////////////////////////////////////////////////////////////////////////
		
		inline void setDT(double dt) 
		{ pdt = dt; }
		
		double getDT(void) const
		{ return pdt; }
		
		////////////////////////////////////////////////////////////////////////
		
		inline double time(void) const
		{ return pTime; }
		
		////////////////////////////////////////////////////////////////////////
		
		uint addComp(steps::sim::CompDef * cdef);
        
		inline uint countComps(void) const
		{ return pComps.size(); }
		
		Comp * comp(uint idx) const;
		
		inline CompPVecCI bgnComp(void) const
		{ return pComps.begin(); }
		inline CompPVecCI endComp(void) const
		{ return pComps.end(); }
		
		////////////////////////////////////////////////////////////////////////
		
		uint addPatch(steps::sim::PatchDef * pdef);
		
		inline uint countPatches(void) const
		{ return pPatches.size(); }
		
		Patch * patch(uint idx) const;
		
		inline PatchPVecCI bgnPatch(void) const
		{ return pPatches.begin(); }
		inline PatchPVecCI endPatch(void) const
		{ return pPatches.end(); }
		
		////////////////////////////////////////////////////////////////////////
		
		Wmrk4 * wmrk4(void) const;
		
		////////////////////////////////////////////////////////////////////////
		
	private:
		
		////////////////////////////////////////////////////////////////////////
		
		steps::sim::StateDef *      pStateDef;
		
		steps::rng::RNG *           pRNG;
		
		
		////////////////////////////////////////////////////////////////////////
		// TIME
		////////////////////////////////////////////////////////////////////////
		
		inline void incTime(double dt)
		{ pTime += dt; }
		
		inline void setTime(double t)
		{ pTime = t; }
		
		inline void resetTime(void)
		{ pTime = 0.0; }
		
		double                      pTime;
		
		/// the maximum time step
		double						pdt;
		
		////////////////////////////////////////////////////////////////////////
		
		Wmrk4 *                     pWmrk4;													
		
		////////////////////////////////////////////////////////////////////////
		// THE MESH ELEMENTS
		////////////////////////////////////////////////////////////////////////
		
		CompPVec                    pComps;
		std::map<steps::sim::CompDef *, Comp *> pCompMap;
		
		PatchPVec                   pPatches;
		
		////////////////////////////////////////////////////////////////////////
		
	};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_WMRK4_SOLVER_CORE_STATE_HPP

// END


