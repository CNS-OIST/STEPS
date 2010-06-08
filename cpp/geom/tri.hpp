////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
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

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_TETMESH_TRI_HPP
#define STEPS_TETMESH_TRI_HPP 1


// STEPS headers.
#include "../common.h"
#include "tet.hpp"
#include "tetmesh.hpp"
#include "tmcomp.hpp"

// STL headers
#include <vector>
#include <ostream>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetmesh)

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
//class TmComp;
class Vertex;
class TmPatch;
class Tet;
class Tetmesh;
//class Tri;

////////////////////////////////////////////////////////////////////////////////

/// A helper class which provides a view on a triangle
/// whose actual data is stored in the TetMesh object.
///
///
/// \warning Methods start with an underscore are not exposed to Python.
class Tri
{

public:

	// TODO: Fix this object before exposure to Python.
	// These objects are not referenced by anything
	// in c++ so ownership is in Python. However, methods like getTet create
	// objects in c++ and cannot be cleaned up by Python.
	// All this data is available in the parent Tetmesh class anyway, so
	// this object may be removed in the future.
	//
	////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

	/// Constructor.
    ///
    /// \param mesh Pointer to the Tetmesh container.
    /// \param tidx Index of the triangle.
    Tri(Tetmesh * mesh, uint tidx);

    /// Destructor.
    ///
    ~Tri(void);

    ////////////////////////////////////////////////////////////////////////
    // TRIANGLE INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns the index of this triangle.
    ///
    /// \return Index of the triangle.
    uint getIdx(void) const
    { return pTidx; }

    /// Returns the area of this triangle.
    ///
    /// \return Area of the triangle.
    double getArea(void) const;

    /// Returns the barycenter of this triangle.
    ///
    /// \return Barycenter of the triangle.
    std::vector<double> getBarycenter(void) const;

    /// Returns the barycenter of this triangle.
    ///
    /// \return Barycenter of the triangle.
    ///
    /// Auxiliary method for internal c++ use.
    double * _getBarycenter(void) const;

    /// Returns the normal of this triangle
    /// by convention points away from the inner tetrahedron.
    ///
    /// \return Normalised triangle.
    std::vector<double> getNorm(void) const;

    /// Returns the normal of this triangle
    /// by convention points away from the inner tetrahedron.
    ///
    /// \return Normalised triangle.
    ///
    /// Auxilliary function to be used internally
    double * _getNorm(void) const;

    /// Returns a pointer to the patch to which this triangle
    /// belongs. Can return 0, if the triangle has not been added to
    /// any compartment.
    ///
    /// \return Pointer to the TmPatch object.
    steps::tetmesh::TmPatch * getPatch(void) const;

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING TETRAHEDRON INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns a Tet object encapsulating 1 of the 2 (possible)
    ///
    /// \param Storing index of the Tet object.
    /// \return Pointer to the Tet object.
    steps::tetmesh::Tet getTet(uint i) const;
    // NOTE: couldn't compile with inline functions
    steps::tetmesh::Tet getTet0(void) const;
    steps::tetmesh::Tet getTet1(void) const;

    /// Returns the inner Tet object.
    ///
    /// \return Pointer to the inner Tet object.
    steps::tetmesh::Tet getInnerTet(void) const;

    /// Returns the outer Tet object.
    ///
    /// \return Pointer to the outer Tet object.
    steps::tetmesh::Tet getOuterTet(void) const;

    ////////////////////////////////////////////////////////////////////////
    /// Returns the index a Tet object encapsulating 1 of the 2 (possible)
    ///
    /// \param Storing index of the Tet object.
    /// \return Index the Tet object.
    int getTetIdx(uint i) const;

    inline int getTet0Idx(void) const
    { return getTetIdx(0); }
    inline int getTet1Idx(void) const
    { return getTetIdx(1); }

    ////////////////////////////////////////////////////////////////////////
    /// Returns the index of the inner Tet object.
    ///
    /// \return Index of the inner Tet object.
    inline int getInnerTetIdx(void) const
    { return getTetIdx(0); }

    /// Returns the index of the outer Tet object.
    ///
    /// \return Index of the outer Tet object.
    inline int getOuterTetIdx(void) const
    { return getTetIdx(1); }

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING VERTEX INFORMATION
    ////////////////////////////////////////////////////////////////////////
    /// Returns the index a vertex encapsulating 1 of the 3 (possible).
    ///
    /// \param Storing index of the vertex.
    /// \return Index the vertex.
    uint getVertexIdx(uint i) const;

    inline uint getVertex0Idx(void) const
    { return getVertexIdx(0); }
    inline uint getVertex1Idx(void) const
    { return getVertexIdx(1); }
    inline uint getVertex2Idx(void) const
    { return getVertexIdx(2); }

private:

    ////////////////////////////////////////////////////////////////////////

	/// Disable the default constructor.
	///
	Tri(void);

	/// A pointer to the parent Tetmesh object
    Tetmesh                   * pTetmesh;

    /// The index of the triangle
    uint                        pTidx;

    /// The 3 vertices of this triangle, by index
    uint                        pVerts[3];

    double 					  * pBaryc;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TRI_HPP

// END
