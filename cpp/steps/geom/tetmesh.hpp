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
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETMESH_TETMESH_HPP
#define STEPS_TETMESH_TETMESH_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STEPS headers.
#include <steps/common.h>
#include <steps/geom/geom.hpp>
#include <steps/geom/tmpatch.hpp>
#include <steps/geom/tmcomp.hpp>

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

/*   The main container class for static tetrahedral meshes.

This class stores the vertices points, tetrahedra and boundary triangles
that comprise the tetrahedral mesh. In addition, it also precomputes
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
    - The total number of tetrahedra in the mesh
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
    - Adding/deleting/moving vertices, triangles and tetrahedra after
      initiation is currently not implemented
*/

class Tetmesh : public steps::wm::Geom
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    Tetmesh(uint nverts, uint ntris, uint ntets);
    Tetmesh(std::vector<double> const & verts, std::vector<uint> const & tets,
    		std::vector<uint> const & tris = std::vector<uint>());
    virtual ~Tetmesh(void);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SETUP
	////////////////////////////////////////////////////////////////////////

    void setVertex(uint vidx, double x, double y, double z);
    void setTri(uint tidx, uint vidx0, uint vidx1, uint vidx2);
    void setTet(uint tidx, uint vidx0, uint vidx1, uint vidx2, uint vidx3);

    void setup(void);
    bool isSetupDone(void) const
    { return pSetupDone; }

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): VERTICES
	////////////////////////////////////////////////////////////////////////

    std::vector<double> getVertex(uint vidx) const;

    inline uint countVertices(void) const
    { return pVertsN; }

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): TRIANGLES
	////////////////////////////////////////////////////////////////////////

    std::vector<uint> getTri(uint tidx) const;

    inline uint countTris(void) const
    { return pTrisN; }

    double getTriArea(uint tidx) const;

    std::vector<double> getTriNorm(uint tidx) const;

    steps::tetmesh::TmPatch * getTriPatch(uint tidx) const;

    void setTriPatch(uint tidx, steps::tetmesh::TmPatch * patch);

    std::vector<int> getTriTetNeighb(uint tidx) const;

    // flip the triangle's inner and outer tetrahedra (not exposed)
    void _flipTriTetNeighb(uint tidx);

    // flip the triangle's vertices and recalculate the normal
    void _flipTriVerts(uint tidx);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): TETRAHEDRA
	////////////////////////////////////////////////////////////////////////

    std::vector<uint>  getTet(uint tidx) const;

    inline uint countTets(void) const
    { return pTetsN; }

    double getTetVol(uint tidx) const;

    steps::tetmesh::TmComp * getTetComp(uint tidx) const;

    void setTetComp(uint tidx, steps::tetmesh::TmComp * comp);

    std::vector<uint> getTetTriNeighb(uint tidx) const;

    std::vector<int> getTetTetNeighb(uint tidx) const;

    // return the index of the tetrahedron that encompasses point;
    // return -1 if point is outside mesh
    // if point is on boundary between two or more tetrahedra,
    // returns first tetrahedron found
    int findTetByPoint(std::vector<double> p) const;

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): MESH
	////////////////////////////////////////////////////////////////////////

    // get the minimal coordinate of the rectangular bounding box
    std::vector<double> getBoundMin(void) const;
    // get the maximal coordinate of the rectangular bounding box
    std::vector<double> getBoundMax(void) const;
    // get the total volume of the mesh
    double getMeshVolume(void) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (C++ INTERNAL)
    ////////////////////////////////////////////////////////////////////////
    
    double * _getVertex(uint vidx) const;
    uint * _getTri(uint tidx) const;
    uint * _getTet(uint tidx) const;
    
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
    /// The triangle normals
    double                            * pTri_norms;
    /// The patch a triangle belongs to
    steps::tetmesh::TmPatch          ** pTri_patches;
    /// The tetrahedra neighbours of each triangle (by index)
    int                               * pTri_tet_neighbours;

    ///////////////////////// DATA: TETRAHEDRA /////////////////////////////
    ///
    /// The total number of tetrahedra in the mesh
    uint                                pTetsN;
    /// The tetrahedra by vertices index
    uint                              * pTets;
    /// The volume of the tetrahedra
    double                            * pTet_vols;
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
