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

#ifndef STEPS_TETMESH_TETMESH_HPP
#define STEPS_TETMESH_TETMESH_HPP 1


// STEPS headers.
#include "../common.h"
#include "geom.hpp"
#include "tmpatch.hpp"
#include "tmcomp.hpp"

// STL headers
#include <vector>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(tetmesh)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class TmPatch;
class TmComp;

////////////////////////////////////////////////////////////////////////////////

template <class T>
bool array_srt_cmp(T ar1[], T ar2[], uint ar_size);

////////////////////////////////////////////////////////////////////////////////

/// The main container class for static tetrahedronl meshes.
/*!
This class stores the vertices points, tetrahedron and boundary triangles
that comprise the tetrahedronl mesh. In addition, it also precomputes
some auxiliary data for the mesh as a whole:

    - Rectangular, axis-aligned bounding box.
    - Overall volume

Auxiliary data is also stored for the tetrahedrons:

    - Volume of each tetrahedron.
    - For each tetrahedron, the indices of its 4 neighbouring tets.
      If there is no neighbour (i.e. if the tetrahedron lies on the
      border), this index will be -1. The sequence of neighbours is
      determined by the following common boundary triangles: (0,1,2);
      (0,1,3); (0,2,3); (1,2,3).
    - For each tetrahedron, the indices of its 4 neighbouring
      boundary triangles. The sequence of neighbours is also
      determined by (0,1,2); (0,1,3); (0,2,3); (1,2,3).
    - The compartment (Comp object) that a tetrahedron belongs to.
      Returns zero pointer if the tetrahedron has not been added to a comp
    - The total number of tetrahedron in the mesh
    - A method of finding which tetrahedron a point given in x,y,z
      coordinates belongs to

And for the triangles:

	- Area of each triangle.
    - Normal vector for each triangle, normalized to length 1.0.
    - For each triangle, the indices of its inside and outside
      neighbouring tetrahedron. If this tetrahedron does not exist
      (because the triangle lies on the outer boundary), this index
      will be -1.
    - The patch (Patch object) that a triangle belongs to. Returns
      zero pointer if triangle has not been added to a patch
    - The total number of triangles in the mesh

And, finally, for the vertices:
    - The total number of vertices in the mesh

NOTES:
    - Adding/deleting/moving vertices, triangles and tetrahedron after
      initiation is currently not implemented

    \warning Methods start with an underscore are not exposed to Python.
*/

class Tetmesh : public steps::wm::Geom
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param nverts Number of vertices.
    /// \param ntets Number of tetrahedrons.
    /// \param ntris Number of triangles
    Tetmesh(uint nverts, uint ntets, uint ntris);

    /// Constructor
    ///
    /// \param verts List of vertices.
    /// \param tets List of tetrahedrons.
    /// \param tris List of triangles.
    Tetmesh(std::vector<double> const & verts, std::vector<uint> const & tets,
    		std::vector<uint> const & tris = std::vector<uint>());

    /// Constructor
    ///
    /// \param verts
    /// TODO: finsih this
    Tetmesh(std::vector<double> const & verts, std::vector<uint> const & tris,
    		std::vector<double> const & tri_areas,
    		std::vector<double> const & tri_norms,
    		std::vector<int> const & tri_tet_neighbs,
    		std::vector<uint> const & tets,
    		std::vector<double> const & tet_vols,
    		std::vector<double> const & tet_barycs,
    		std::vector<uint> const & tet_tri_neighbs,
    		std::vector<int> const & tet_tet_neighbs);

    /// Destructor
    virtual ~Tetmesh(void);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SETUP
	////////////////////////////////////////////////////////////////////////
    /// Setup a vertex.
    ///
    /// \param vidx Index of the vertex.
    /// \param x Coordinate x.
    /// \param y Coordinate y.
    /// \param z Coordinate z.
    void setVertex(uint vidx, double x, double y, double z);

    /// Setup a triangle.
    ///
    /// \param tidx Index of the triangle.
    /// \param vidx0 Index of vertex 0 that forms the triangle.
    /// \param vidx1 Index of vertex 1 that forms the triangle.
    /// \param vidx2 Index of vertex 2 that forms the triangle.
    void setTri(uint tidx, uint vidx0, uint vidx1, uint vidx2);

    /// Setup a tetrahedron.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param vidx0 Index of vertex 0 that forms the tetrahedron.
    /// \param vidx1 Index of vertex 1 that forms the tetrahedron.
    /// \param vidx2 Index of vertex 2 that forms the tetrahedron.
    /// \param vidx3 Index of vertex 3 that forms the tetrahedron.
    void setTet(uint tidx, uint vidx0, uint vidx1, uint vidx2, uint vidx3);

    /// Setup the Temesh.
    void setup(void);

    /// Check if the Temesh is set up.
    ///
    /// \return True if the Temesh is set up;
    ///         False if else.
    bool isSetupDone(void) const
    { return pSetupDone; }

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): VERTICES
	////////////////////////////////////////////////////////////////////////

    /// Return the coordinates of a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    std::vector<double> getVertex(uint vidx) const;

    /// Count the vertices in the Temesh.
    ///
    /// \return Number of the vertices.
    inline uint countVertices(void) const
    { return pVertsN; }

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): TRIANGLES
	////////////////////////////////////////////////////////////////////////

    /// Return the triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Indices of the vertices form the triangle.
    std::vector<uint> getTri(uint tidx) const;

    /// Count the triangles in the Temesh.
    ///
    /// \return Number of the triangles.
    inline uint countTris(void) const
    { return pTrisN; }

    /// Return the area of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Area of the triangle.
    double getTriArea(uint tidx) const;

    /// Return the barycenter of triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Barycenter of the triangle.
    std::vector<double> getTriBarycenter(uint tidx) const;


    /// Return the normalised triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Coordinate of the normalised vertices form the triangle.
    std::vector<double> getTriNorm(uint tidx) const;

    //getTriBarycenter
    /// Return the patch which a triangle associated to.
    ///
    /// \param tidx Index of the triangle.
    /// \return Pointer to the patch.

    steps::tetmesh::TmPatch * getTriPatch(uint tidx) const;
    ///Set the patch which a triangle belongs to.
    ///
    /// \param tidx Index of the triangle.
    /// \param patch Pointer to the associated patch.

    void setTriPatch(uint tidx, steps::tetmesh::TmPatch * patch);
    ///Return the tetrahedron neighbors of a triangle by its index.
    ///
    /// \param tidx Index of the triangle.
    /// \return Vector of the tetrahedron neighbors.

    std::vector<int> getTriTetNeighb(uint tidx) const;
    
    /// Return the triangles which form the surface boundary of the mesh.
    /// \return Vector of the triangle boundary.
    // Weiliang 2010.02.02
    std::vector<int> getTriBoundary(void) const;

    /// Flip the triangle's inner and outer tetrahedron.
    ///
    /// \param tidx Index of the triangle.
    void _flipTriTetNeighb(uint tidx);

    /// Flip the triangle's vertices and recalculate the normal.
    ///
    /// \param Index of the triangle.
    void _flipTriVerts(uint tidx);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): TETRAHEDRA
	////////////////////////////////////////////////////////////////////////
    /// Return a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the indices of triangles which form the tetrahedron.
    std::vector<uint>  getTet(uint tidx) const;
    /// Count the number of tetrahedrons.
    ///
    /// \return Number of tetrahedrons.
    inline uint countTets(void) const
    { return pTetsN; }
    /// Return the volume of a tetrahedron by its index.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Volume of the tetrahedron.

    double getTetVol(uint tidx) const;

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
    /// \return Quality RER of the tetrahedron.
    double getTetQualityRER(uint tidx) const;

    /// Return the barycentre of the tetrahedron in x,y,z coordinates
    /// \param tidx Index of the tetrahedron
    /// \return Barycentre of the tetrahedron
	std::vector<double> getTetBarycenter(uint tidx) const;

    /// Return the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Pointer to the compartment object.

    steps::tetmesh::TmComp * getTetComp(uint tidx) const;
    ///Set the compartment which a tetrahedron with index tidx belongs to.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \param comp Pointer to the compartment object.

    void setTetComp(uint tidx, steps::tetmesh::TmComp * comp);
    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the triangle neighbors.

    std::vector<uint> getTetTriNeighb(uint tidx) const;
    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Vector of the tetrahedron neighbors.

    std::vector<int> getTetTetNeighb(uint tidx) const;

    /// Find a tetrahedron which encompasses a given point.
    /// Return the index of the tetrahedron that encompasses point;
    ///  return -1 if point is outside mesh;
    /// if point is on boundary between two or more tetrahedron,
    /// returns first tetrahedron found.
    /// \param p A point given by its coordinates.
    /// \return ID of the found tetrahedron.

    int findTetByPoint(std::vector<double> p) const;

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): MESH
	////////////////////////////////////////////////////////////////////////

    /// Return the minimal coordinate of the rectangular bounding box.
    ///
    /// \return Minimal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMin(void) const;
    /// Return the maximal coordinate of the rectangular bounding box.
    ///
    /// \return Maximal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMax(void) const;
    /// Return the total volume of the mesh.
    ///
    /// \return Volume of the mesh.
    double getMeshVolume(void) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (C++ INTERNAL)
    ////////////////////////////////////////////////////////////////////////

    /// Return a vertex with index vidx.
    ///
    /// \param vidx Index of the vertex.
    /// \return Coordinates of the vertex.
    double * _getVertex(uint vidx) const;

    /// Return a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return List of the vertices form the triangle.
    uint * _getTri(uint tidx) const;

    /// Return a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return List of the vertices form the tetrahedron.
    uint * _getTet(uint tidx) const;

    ///Return the tetrahedron neighbors of a triangle with index tidx.
    ///
    /// \param tidx Index of the triangle.
    /// \return Array of the tetrahedron neighbors.
    int * _getTriTetNeighb(uint tidx) const;

    ///Return the triangle neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the triangle neighbors.
    uint * _getTetTriNeighb(uint tidx) const;

    ///Return the tetrahedron neighbors of a tetrahedron with index tidx.
    ///
    /// \param tidx Index of the tetrahedron.
    /// \return Array of the tetrahedron neighbors.
    int * _getTetTetNeighb(uint tidx) const;

    /// Return the normalised triangle with index tidx
    ///
    /// \param tidx Index of the triangle.
    /// \return Array of Coordinate of the normalised vertices form the triangle.
    double * _getTriNorm(uint tidx) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    bool                                pSetupDone;

    ///////////////////////// DATA: VERTICES ///////////////////////////////
    ///
    /// The total number of vertices in the mesh
    uint                                pVertsN;
    /// The vertices by x,y,z coordinates
    double                            * pVerts;

    ///////////////////////// DATA: TRIANGLES //////////////////////////////
    ///
    /// The total number of triangles in the mesh
    uint                                pTrisN;
    /// The triangles by vertices index
    uint                              * pTris;
    /// Array available for user-supplied triangle data for 2nd constructor
    uint                              * pTris_user;
    /// The areas of the triangles
    double                            * pTri_areas;
    /// The triangle barycenters
    double                            * pTri_barycs;
    /// The triangle normals
    double                            * pTri_norms;
    /// The patch a triangle belongs to
    steps::tetmesh::TmPatch          ** pTri_patches;
    /// The tetrahedron neighbours of each triangle (by index)
    int                               * pTri_tet_neighbours;

    ///////////////////////// DATA: TETRAHEDRA /////////////////////////////
    ///
    /// The total number of tetrahedron in the mesh
    uint                                pTetsN;
    /// The tetrahedron by vertices index
    uint                              * pTets;
    /// The volume of the tetrahedron
    double                            * pTet_vols;
	/// The barycentres of the tetrahedra
	double                            * pTet_barycentres;
    /// The compartment a tetrahedron belongs to
    steps::tetmesh::TmComp           ** pTet_comps;
    /// The triangle neighbours of each tetrahedron (by index)
    uint                              * pTet_tri_neighbours;
    /// The tetrahedron neighbours of each tetrahedron (by index)
    int                               * pTet_tet_neighbours;



    ////////////////////////////////////////////////////////////////////////

    /// Information about the minimal and maximal boundary values
    ///
    double                      pXmin;
    double                      pXmax;
    double                      pYmin;
    double                      pYmax;
    double                      pZmin;
    double                      pZmax;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(tetmesh)
END_NAMESPACE(steps)

#endif
// STEPS_TETMESH_TETMESH_HPP

// END
