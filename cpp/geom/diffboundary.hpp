////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2011ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETMESH_DIFFBOUNDARY_HPP
#define STEPS_TETMESH_DIFFBOUNDARY_HPP 1


// STEPS headers.
#include "../common.h"
#include "diffboundary.hpp"
#include "tetmesh.hpp"
#include "../model/surfsys.hpp"
#include "tmcomp.hpp"
#include "geom.hpp"

// STL headers
#include <vector>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetmesh)

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;
class DiffBoundary;

// Auxiliary declarations.
typedef DiffBoundary *                          DiffBoundaryP;
typedef std::map<std::string, DiffBoundaryP>    DiffBoundaryPMap;
typedef DiffBoundaryPMap::iterator              DiffBoundaryPMapI;
typedef DiffBoundaryPMap::const_iterator        DiffBoundaryPMapCI;

typedef std::vector<DiffBoundaryP>              DiffBoundaryPVec;
typedef DiffBoundaryPVec::iterator              DiffBoundaryPVecI;
typedef DiffBoundaryPVec::const_iterator        DiffBoundaryPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of diffusion boundary triangles of a Tetmesh.
///
/// \warning Methods start with an underscore are not exposed to Python.
class DiffBoundary
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the DiffBoundary.
    /// \param container Pointer to the Tetmesh container.
    /// \param tris A sequence of triangles (by index) as a vector
	///             of unsigned integers which is represented as
    ///             a sequence of positive integer values) in Python.

    ///
    /// This is the constructor for the tetmesh (tetrahedron mesh) namespace.
    DiffBoundary(std::string const & id, Tetmesh * container,
       	 std::vector<uint> const & tris);

    /// Destructor.
    virtual ~DiffBoundary(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

	/// Return the patch id.
    ///
    /// \return ID of the diffusion boundary.
	std::string getID(void) const
	{ return pID; }

	/// Set or change the diffusion boundary id.
    ///
    /// \param id ID of the diffusion boundary.
	void setID(std::string const & id);

	/// Return a pointer to the geometry container object.
    ///
    /// \return Pointer to the parent geometry container.
	steps::tetmesh::Tetmesh * getContainer(void) const
	{ return pTetmesh; }

    /// Return whether triangles (specified by index) are inside this patch.
    ///
    /// \param tri List of indices of triangles.
    /// \return Results of whether the triangles are inside the patch.
    std::vector<bool> isTriInside(std::vector<uint> tri) const;

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<uint> getAllTriIndices(void) const
    { return pTri_indices; }

	/// Return the compartments this diffusion boundary connects
    ///
    /// \return List of the two compartments.
    inline std::vector<steps::wm::Comp *>  getComps(void) const
	{ return pComps; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<uint> const & _getAllTriIndices(void) const
    { return pTri_indices; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
	/// Called if Python object deleted, or from del method in parent object.
	/// Will only be called once
	void _handleSelfDelete(void);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;

	std::vector<steps::wm::Comp *>      pComps;

	steps::tetmesh::Tetmesh 		  * pTetmesh;
	std::vector<uint>                   pTri_indices;
	uint								pTrisN;



////////////////////////////////////////////////////////////////////////////////

};

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_DIFFBOUNDARY_HPP

// END
