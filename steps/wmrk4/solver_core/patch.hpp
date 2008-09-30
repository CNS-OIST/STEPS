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

#ifndef STEPS_WMRK4_SOLVER_CORE_PATCH_HPP
#define STEPS_WMRK4_SOLVER_CORE_PATCH_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include <steps/common.h>
#include <steps/sim/shared/patchdef.hpp>

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Patch;
class Comp;

// Auxiliary declarations.
typedef Patch *                         PatchP;
typedef std::vector<PatchP>             PatchPVec;
typedef PatchPVec::iterator             PatchPVecI;
typedef PatchPVec::const_iterator       PatchPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Patch
	{
		
	public:
		
		////////////////////////////////////////////////////////////////////////
		// OBJECT CONSTRUCTION & DESTRUCTION
		////////////////////////////////////////////////////////////////////////
		
		Patch(steps::sim::PatchDef * patchdef, Comp * icomp, Comp * ocomp);
		~Patch(void);
		
		void reset(void);
		
		////////////////////////////////////////////////////////////////////////
		// DATA ACCESS
		////////////////////////////////////////////////////////////////////////
		
		inline steps::sim::PatchDef * def(void) const
		{ return pPatchDef; }
		
		inline double area(void) const
		{ return pPatchDef->area(); }
		
		////////////////////////////////////////////////////////////////////////
		
		inline double * pools(void) const													
		{ return pPoolCount; }
		
		inline uint * flags(void) const
		{ return pPoolFlags; }
		
		inline uint * srflags(void) const
		{ return pSReacFlags; }
		
		static const uint CLAMPED = 1;
		
		inline bool clamped(uint lidx) const
		{ return pPoolFlags[lidx] & CLAMPED; }
		void setClamped(uint lidx, bool clamp);
		
		static const uint INACTIVATED = 1;
		
		inline bool active(uint lidx) const
		{ return !(pSReacFlags[lidx] & INACTIVATED); }
		void setActive(uint lidx, bool active);
		
		////////////////////////////////////////////////////////////////////////
		
		inline Comp * iComp(void) const
		{ return pIComp; }
		
		inline Comp * oComp(void) const
		{ return pOComp; }
		
		////////////////////////////////////////////////////////////////////////
		
	private:
		
		////////////////////////////////////////////////////////////////////////
		
		steps::sim::PatchDef *      pPatchDef;
		
		/// Numbers of molecules -- stored as doubles for numerical solver.					
		double *                    pPoolCount;															
		/// Flags on these pools -- stored as machine word flags.
		uint *                      pPoolFlags;
		/// Flags on sreactions, stored as machine word srflags
		uint *						pSReacFlags;															
		
		Comp *                      pIComp;
		Comp *                      pOComp;
		
		////////////////////////////////////////////////////////////////////////
		
		// Disable constructors.
		Patch(void);
		Patch(Patch const &);
		
		////////////////////////////////////////////////////////////////////////
		
	};

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_WMRK4_SOLVER_CORE_PATCH_HPP

// END


