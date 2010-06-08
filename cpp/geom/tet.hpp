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

#ifndef STEPS_TETMESH_TET_HPP
#define STEPS_TETMESH_TET_HPP 1


// STEPS headers.
#include "../common.h"
#include "tri.hpp"
#include "tetmesh.hpp"
#include "tmcomp.hpp"
#include "../rng/rng.hpp"


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
class Tri;

////////////////////////////////////////////////////////////////////////////////

/// Helper class which provides a view on a tetrahedron
/// whose actual data is stored in the TetMesh object.
///
/// \todo: implement getQualityAR
/// \warning Methods start with an underscore are not exposed to Python.
class Tet
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
    /// \param mesh The Tetmesh container.
    /// \param tidx Index of the tetrahedron.
    Tet(Tetmesh * mesh, uint tidx);

    /// Destructor.
    ~Tet(void);

    ////////////////////////////////////////////////////////////////////////
    // TETRAHEDRON INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns the index of this tetrahedron.
    ///
    /// \return Index of the tetrahedron.
    uint getIdx(void) const
    { return pTidx; }

    /// Returns the volume of this tetrahedron.
    ///
    /// \return Volume of the tetrahedron.
    double getVol(void) const;

    /// Returns the barycenter of this tetrahedron.
    ///
    /// \return Coordinates of the barycenter.
    std::vector<double> getBarycenter(void) const;

    /// Returns the barycenter of this tetrahedron.
    ///
    /// \return Coordinates of the barycenter.
    ///
    /// Auxiliary method for internal c++ use.
    double * _getBarycenter(void) const;

    // Computes the quality of the tetrahedron.
    //
    // This method uses the aspect ratio (AR) metric for tetrahedron
    // quality, given by dividing the length of the longest edge with
    // the smallest altitude. The smaller this value, the more regular
    // the tetrahedron.
    //


    //double getQualityAR(void) const;

    /// Computes the quality of the tetrahedron.
    ///
    /// This method uses the radius-edge ratio (RER) metric for tetrahedron
    /// quality, given by dividing the radius of the tetrahedron's
    /// circumsphere with the length of the shortest edge.
    ///
    /// The smaller this value, the more regular the tetrahedron. The
    /// lowest possible value of this metric is given by computing the
    /// RER for a fully regular tetrahedron:
    ///
    ///    Q = sqrt(6)/4 ~ 0.612
    ///
    /// This is a slightly weaker metric than getQualityAR, because
    /// certain slivers (degenerate tetrahedrons) can still have a fairly
    /// small value.
    ///
    /// \return Quality of the tetrahedron.
    double getQualityRER(void) const;

    ////////////////////////////////////////////////////////////////////////

    /// Returns a pointer to the compartment to which this tetrahedron
    /// belongs.
    ///
    /// Can return 0, if the tetrahedron has not been added to
    /// any compartment.
    /// \return Pointer to the compartment.
    steps::tetmesh::TmComp * getComp(void) const;

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING TETRAHEDRON INFORMATION
    ////////////////////////////////////////////////////////////////////////

    /// Returns a Tet object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tet object.
    /// \return Pointer to the Tet object.
    Tet getTet(uint i) const;

    /// Return the Tet object 0.
    ///
    /// \return Pointer to the Tet object.
    inline Tet getTet0(void) const
    { return getTet(0); }

    /// Return the Tet object 1.
    ///
    /// \return Pointer to the Tet object.
    inline Tet getTet1(void) const
    { return getTet(1); }

    /// Return the Tet object 2.
    ///
    /// \return Pointer to the Tet object.
    inline Tet getTet2(void) const
    { return getTet(2); }

    /// Return the Tet object 3.
    ///
    /// \return Pointer to the Tet object.
    inline Tet getTet3(void) const
    { return getTet(3); }

    ////////////////////////////////////////////////////////////////////////
    /// Returns the index a Tet object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tet object.
    /// \return Index of the Tet object.
    int getTetIdx(uint i) const;

    /// Return the Tet object 0.
    ///
    /// \return Index of the Tet object.
    inline int getTet0Idx(void) const
    { return getTetIdx(0); }

    /// Return the Tet object 1.
    ///
    /// \return Index of the Tet object.
    inline int getTet1Idx(void) const
    { return getTetIdx(1); }

    /// Return the Tet object 2.
    ///
    /// \return Index of the Tet object.
    inline int getTet2Idx(void) const
    { return getTetIdx(2); }

    /// Return the Tet object 3.
    ///
    /// \return Index of the Tet object.
    inline int getTet3Idx(void) const
    { return getTetIdx(3); }

    ////////////////////////////////////////////////////////////////////////
    /// Returns the distance between the barycenter
    /// and a Tet object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tet object.
    /// \return Distance between the barycenter and the Tet object.
    double getTetDist(uint i) const;

    /// Returns the distance between the barycenter and the Tet object 0.
    inline double getTet0Dist(void) const
    { return getTetDist(0); }

    /// Returns the distance between the barycenter and the Tet object 1.
    inline double getTet1Dist(void) const
    { return getTetDist(1); }

    /// Returns the distance between the barycenter and the Tet object 2.
    inline double getTet2Dist(void) const
    { return getTetDist(2); }

    /// Returns the distance between the barycenter and the Tet object 3.
    inline double getTet3Dist(void) const
    { return getTetDist(3); }

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING TRIANGLE INFORMATION
    ////////////////////////////////////////////////////////////////////////
    /// Returns a Tri object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tri object.
    /// \return Pointer to the Tri object.
    Tri getTri(uint i) const;
    // NOTE: couldn't compile with inline functions
    Tri getTri0(void) const;
    Tri getTri1(void) const;
    Tri getTri2(void) const;
    Tri getTri3(void) const;

    ////////////////////////////////////////////////////////////////////////
    /// Returns the index of a Tri object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tri object.
    /// \return Index of the Tri object.
    uint getTriIdx(uint i) const;

    inline uint getTri0Idx(void) const
    { return getTriIdx(0); }
    inline uint getTri1Idx(void) const
    { return getTriIdx(1); }
    inline uint getTri2Idx(void) const
    { return getTriIdx(2); }
    inline uint getTri3Idx(void) const
    { return getTriIdx(3); }

    /// Returns the distance between the barycenter
    /// and a Tri object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tri object.
    /// \return Distance between the barycenter and the Tri object.
    double getTriDist(uint i) const;

    inline double getTri0Dist(void) const
    { return getTriDist(0); }
    inline double getTri1Dist(void) const
    { return getTriDist(1); }
    inline double getTri2Dist(void) const
    { return getTriDist(2); }
    inline double getTri3Dist(void) const
    { return getTriDist(3); }

    ////////////////////////////////////////////////////////////////////////
    /// Returns the area of a Tri object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Tri object.
    /// \return Area of the Tri object.
    double getTriArea(uint i) const;

    inline double getTri0Area(void) const
    { return getTriArea(0); }
    inline double getTri1Area(void) const
    { return getTriArea(1); }
    inline double getTri2Area(void) const
    { return getTriArea(2); }
    inline double getTri3Area(void) const
    { return getTriArea(3); }

    ////////////////////////////////////////////////////////////////////////
    // NEIGHBOURING VERTEX INFORMATION
    ////////////////////////////////////////////////////////////////////////
    /// Returns the index of a Vertex object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Vertex object.
    /// \return Index of the Vertex.
    uint getVertexIdx(uint i) const;

    inline uint getVertex0Idx(void) const
    { return getVertexIdx(0); }
    inline uint getVertex1Idx(void) const
    { return getVertexIdx(1); }
    inline uint getVertex2Idx(void) const
    { return getVertexIdx(2); }
    inline uint getVertex3Idx(void) const
    { return getVertexIdx(3); }

    /// Returns the coordinates of a Vertex object encapsulating 1 of the 4 (possible).
    ///
    /// \param i Encapuslating index of the Vertex object.
    /// \return Coordinates of the Vertex.
    std::vector<double> getVertex(uint i) const;

    /// Returns the coordinates of a Vertex object encapsulating 1 of the 4 (possible).
    ///
    /// Auxiliary function for internal use.
    /// \param i Encapuslating index of the Vertex object.
    /// \return Coordinates of the Vertex.
    double * _getVertex(uint i) const;

    inline uint getVertex0(void) const
    { return getVertexIdx(0); }
    inline uint getVertex1(void) const
    { return getVertexIdx(1); }
    inline uint getVertex2(void) const
    { return getVertexIdx(2); }
    inline uint getVertex3(void) const
    { return getVertexIdx(3); }


    ////////////////////////////////////////////////////////////////////////
    // OTHER FUNCTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Returns true if the point is inside the tetrahedron.
    ///
    /// \param p coordinates of a point.
    /// \return True if the point is inside the tetrahedron;
    ///         False if else.
    bool isInside(std::vector<double> p) const;

    /// Generate a number of random points in this tetrahedron. The
    /// default number is 1. Returns a pointer to an array of size N * 3
    /// of doubles.
    ///
    /// \param r Random number generator.
    /// \param n Number of points to be generated.
    /// \return Coordinates of generated points.
    std::vector<double> getRanPnt(steps::rng::RNG * r, uint n = 1) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    /// Disable the default constructor.
    ///
    Tet(void);

    /// A pointer to the parent Tetmesh object
    Tetmesh                   * pTetmesh;

    /// The index of this tetrahedron
    uint                        pTidx;

    /// The four vertices of this tetrahedron, by index
    uint                        pVerts[4];

    // The Barycentre- currently calculated each time _getBarycentre is called
    double                     * pBaryc;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TET_HPP

// END
