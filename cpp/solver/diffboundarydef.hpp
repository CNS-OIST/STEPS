////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011�Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006�University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS�is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS�is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.�If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////


#ifndef STEPS_SOLVER_DIFFBOUNDARYDEF_HPP
#define STEPS_SOLVER_DIFFBOUNDARYDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "../common.h"
#include "statedef.hpp"
#include "../geom/diffboundary.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(solver)

// Forwards declarations
class DiffBoundarydef;

// Auxiliary declarations.
typedef DiffBoundarydef *                      DiffBoundaryDefP;
typedef std::vector<DiffBoundaryDefP>          DiffBoundaryDefPVec;
typedef DiffBoundaryDefPVec::iterator          DiffBoundaryDefPVecI;
typedef DiffBoundaryDefPVec::const_iterator    DiffBoundaryDefPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Defined diffusion boundary object.
class DiffBoundarydef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param d Pointer to the associated Diff boundary object.
	DiffBoundarydef(Statedef * sd, uint idx, steps::tetmesh::DiffBoundary * db);

    /// Destructor
	~DiffBoundarydef(void);

	void setup(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION BOUNDARY
    ////////////////////////////////////////////////////////////////////////

	/// Return the global index of this diffusion boundary.
	inline uint gidx(void) const
	{ return pIdx; }

    /// Return the name of this diffusion boundary.
	std::string const name(void) const;

	inline std::vector<uint> tris(void) const
	{ return pTris; }

	inline uint compa(void) const
	{ return pCompA; }
	inline uint compb(void) const
	{ return pCompB; }


	////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

	Statedef                          * pStatedef;

	bool 								pSetupdone;

	// The global index of this diffusion boundary
	uint                                pIdx;

	// The string identifier of this diffusion rule
	std::string                         pName;

	// List of all the triangles

	std::vector<uint> 					pTris;

	// Diffboundarydef will have a setup() to copy the Compdef pointers
	//Lets store pointers
	uint 								pCompA;
	uint					 			pCompB;

	// The pointer to the well-mixed comps is stored, but should not be used
	// only here so it's available during setup.
	steps::wm::Comp * 					pCompA_temp;
	steps::wm::Comp * 					pCompB_temp;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(solver)
END_NAMESPACE(steps)

#endif
// STEPS_SOLVER_DIFFBOUNDARYDEF_HPP

// END
