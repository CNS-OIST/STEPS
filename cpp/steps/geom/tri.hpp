////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
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

#ifndef STEPS_TETMESH_TRI_HPP
#define STEPS_TETMESH_TRI_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/geom/tet.hpp>
#include <steps/geom/tetmesh.hpp>
#include <steps/geom/tmcomp.hpp>

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

/// Class Tri is a helper class which provides a view on a triangle
/// whose actual data is stored in the TetMesh object.
///
class Tri
{

public:

	////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

	/// Constructor.
    ///
    Tri(Tetmesh * mesh, uint tidx);

    /// Destructor.
    ///
    ~Tri(void);

    ////////////////////////////////////////////////////////////////////////
    // TRIANGLE INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns the index of this triangle.
    ///
    uint getIdx(void) const
    { return pTidx; }

    /// Returns the area of this triangle.
    ///
    double getArea(void) const;

    /// Returns the barycenter of this triangle.
    ///
    std::vector<double> getBarycenter(void) const;

    /// Returns the normal of this triangle
    /// by convention points away from the inner tetrahedron
    std::vector<double> getNorm(void) const;

    /// Returns a pointer to the patch to which this triangle
    /// belongs. Can return 0, if the triangle has not been added to
    /// any compartment
    ///
    steps::tetmesh::TmPatch * getPatch(void) const;

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING TETRAHEDRON INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns a Tet object encapsulating 1 of the 2 (possible)
    steps::tetmesh::Tet getTet(uint i) const;
    // NOTE: couldn't compile with inline functions
    steps::tetmesh::Tet getTet0(void) const;
    steps::tetmesh::Tet getTet1(void) const;

    steps::tetmesh::Tet getInnerTet(void) const;
    steps::tetmesh::Tet getOuterTet(void) const;

    ////////////////////////////////////////////////////////////////////////

    int getTetIdx(uint i) const;

    inline int getTet0Idx(void) const
    { return getTetIdx(0); }
    inline int getTet1Idx(void) const
    { return getTetIdx(1); }

    ////////////////////////////////////////////////////////////////////////

    inline int getInnerTetIdx(void) const
    { return getTetIdx(0); }
    inline int getOuterTetIdx(void) const
    { return getTetIdx(1); }

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING VERTEX INFORMATION
    ////////////////////////////////////////////////////////////////////////

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

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TRI_HPP

// END
