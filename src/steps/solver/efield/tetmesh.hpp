/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */


#ifndef STEPS_SOLVER_EFIELD_TETMESH_HPP
#define STEPS_SOLVER_EFIELD_TETMESH_HPP 1

// STL headers.
#include <iostream>
#include <map>
#include <set>
#include <string>
#include <vector>
#include <queue>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/efield/vertexconnection.hpp"
#include "steps/solver/efield/vertexelement.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps{
 namespace solver {
 namespace efield {
using namespace std;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class TetMesh;

////////////////////////////////////////////////////////////////////////////////

/// Objects of this class are used to store tetrahedrons in a format
/// that makes it easy to find whether a certain quadruplet of vertices
/// corresponds to a tetrahedron.
///
/// \author Stefan Wils
///
class TetStub
{

public:

    /// Constructor.
    ///
    TetStub(uint v1, uint v2, uint v3, uint v4);

    /// Constructor.
    ///
    TetStub(uint * v);

    /// Comparison operator -- required for storing objects of this
    /// class in a set. Implements strict weak ordering.
    ///
    bool operator< (TetStub const & t) const;

    /// Output stream operator.
    ///
    // friend std::ostream & operator<< (std::ostream & os, steps::solver::efield::TetStub const &);



private:

    /// Vertices on the edges of the tetrahedron -- sorted from
    /// small to large.
    uint        pSortedVerts[4];

};

////////////////////////////////////////////////////////////////////////////////

// I.H. moved from .cpp file
struct ConnStub
{

    ConnStub(VertexElement * vertex1, VertexElement * vertex2)
    : fVertex1(vertex1)
    , fVertex2(vertex2)
    {
        if (fVertex1 < fVertex2)
        {
            VertexElement * tmp = fVertex1;
            fVertex1 = fVertex2;
            fVertex2 = tmp;
        }
    }

    bool operator< (ConnStub const & c) const
    {
        if (fVertex1 < c.fVertex1) return true;
        if (fVertex1 > c.fVertex1) return false;
        if (fVertex2 < c.fVertex2) return true;
        return false;
    }
    /*
    // I.H. 24/11/09  Isn't this necessary for a set??
    bool operator == (ConnStub const & c) const
    {
        if (fVertex1 != c.fVertex1) return false;
        if (fVertex2 != c.fVertex2) return false;
        return true;
    }
    */
    VertexElement * fVertex1;
    VertexElement * fVertex2;

};

////////////////////////////////////////////////////////////////////////////////

// Auxiliary declarations.
typedef std::set<TetStub>               TetStubSet;
typedef TetStubSet::iterator            TetStubSetI;
typedef TetStubSet::const_iterator      TetStubSetCI;

////////////////////////////////////////////////////////////////////////////////

/// Class TetMesh stores a copy of the mesh in a format which is suitable
/// for the e-field calculations. This means it stores vertices as separate
/// objects (of class VertexElement) and also explicitly stores the edges
/// (in objects of class VertexConnection).
///
/// This class used to be derived from Mesh (now deleted -- should be
/// available in the history of the Subversion archive). We've merged
/// them because there is no reason for more than one mesh type in the
/// near future and because it was a mess anyway.
///
/// \author Robert Cannon
///
class TetMesh
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    TetMesh(uint nv, double * vpos,
            uint ntr, uint * trivi,
            uint ntet, uint * tetvi);
    ~TetMesh(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    /// Called by the EField constructor after all the triangles and
    /// tetrahedrons have been specified. It extracts all unique
    /// vertex-vertex connections by looping over all tetrahedrons.
    ///
    /// Originally from TetMesh.
    ///
    void extractConnections(void);

    /// Called by the EField constructor after extractConnections().
    /// This method basically loops over all triangles in the mesh,
    /// computes their area, and adds one third of this area to each
    /// vertex of that triangle. In the end, the total surface area
    /// associated with each vertex therefore contains contributions
    /// from all triangles that share this vertex.
    ///
    /// Originally from TetMesh.
    ///
    void allocateSurface(void);

    ////////////////////////////////////////////////////////////////////////

    /// This method (re-?)assigns all VertexElements new indices, based
    /// on their current order in the vector stored here. This order
    /// might change.
    ///
    /// Originally from Mesh.
    ///
    void reindexElements(void);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VERTICES
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of vertices in the mesh.
    ///
    /// Originally from Mesh.
    ///
    inline uint countVertices(void) const
    { return pElements.size(); }

    /// Originally from Mesh.
    ///
    inline VertexElement * getVertex(uint i) const
    { return pElements[i]; }

    /// For a given vertex, it constructs a list of tetrahedrons that
    /// include this vertex. The tetrahedra are returned in a special
    /// format:
    ///
    /// <UL>
    /// <LI> For each tetrahedron, only the three 'other' points are
    ///      retured. In other words, the vertex for which this method
    ///      is called is itself not included. So each array stored in
    ///      the vector has size 3.
    /// <LI> For each tetrahedron, these three 'other' points are not
    ///      indices in VertexElement's list of neighbours.
    /// </UL>
    ///
    /// This information is used by class TetCoupler to compute coupling
    /// constants.
    ///
    /// Originally from TetMesh.
    ///
    std::vector<std::vector<uint> >  getNeighboringTetrahedra(VertexElement *);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: TRIANGLES
    ////////////////////////////////////////////////////////////////////////

    /// Originally from TetMesh.
    ///
    inline uint getNTri(void) const
    { return pNTri; }

    /// Originally from TetMesh.
    ///
    inline uint * getTriangle(uint i) const
    { return pTriangles + (3 * i); }

    /// Originally from TetMesh.
    ///
    inline uint getTriangleVertex(uint itr, uint iv) const
    { return pTriangles[(3 * itr) + iv]; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: TETRAHEDRONS
    ////////////////////////////////////////////////////////////////////////

    inline uint getNTet(void) const
    { return pNTet; }

    inline uint * getTetrahedron(uint i) const
    { return pTetrahedrons + (4 * i); }

    /// Originally from TetMesh.
    ///
    inline uint getTetrahedronVertex(uint itr, uint iv) const
    { return pTetrahedrons[(4 * itr) + iv]; }

    ////////////////////////////////////////////////////////////////////////

    /// Originally from Mesh.
    ///
    void applySurfaceCapacitance(double);

    void applyTriCapacitance(uint tidx, double cm);


    /// Originally from Mesh.
    ///
    void applyConductance(double);

    /// Originally from Mesh.
    ///
    double getTotalCapacitance(void);

    /// Originally from Mesh.
    ///
    double getTotalArea(void);

    /// Originally from Mesh.
    /// Iain: big changes here
    ///
    void axisOrderElements(uint opt_method, std::string const & opt_file_name ="", double search_percent=100.0);

    void saveOptimal(std::string const & opt_file_name);

    void fill_ve_vec(set<VertexElement*> & veset, vector<VertexElement*> & vevec, queue<VertexElement*> & vequeue, uint ncons, VertexElement ** nbrs);

    /// Originally from Mesh.
    ///
    std::vector<double> centerOfMass(void);

    /// Originally from Mesh.
    /// Iain: defunct: replaced by 'walking' method in fill_ve_vec function
    ///
    void mainEvec(uint, double**, uint*, double*, double*);

    /// Originally from Mesh.
    ///
    uint ncon(void)
    { return pConnections.size(); }

    /// Originally from Mesh.
    ///
    VertexConnection* getConnection(uint i)
    { return pConnections[i]; }

    std::vector<uint> getVertexPermutation(void);

    ////////////////////////////////////////////////////////////////////////
    // FROM TETMESH
    ////////////////////////////////////////////////////////////////////////

    /// Originally from TetMesh.
    ///
    void reconstruct(void);

    /// Originally from TetMesh.
    ///
    void reordered(void);

    // Just temporary functions to have a look at the matrix
    //void displayMatrix(uint);
    //void savematrix(void);

private:

    ////////////////////////////////////////////////////////////////////////
    // AUXILIARY OBJECT CONSTRUCTION & SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Adds a (unique) connection between two vertex elements. This
    /// is called from TetMesh::extractConnections.
    ///
    /// swils -- this will currently lead to a lot of vector resizing
    /// and subsequent copying going on during initialization. If this
    /// would become a problem, maybe look for a way to 'predict' the
    /// number of connections and then allocate a suitably sized vector
    /// in advance?
    ///
    /// Originally from Mesh.
    ///
    VertexConnection * newConnection(VertexElement *, VertexElement *);

    ////////////////////////////////////////////////////////////////////////

    ////////////////////////////////////////////////////////////////////////
    // COPIED FROM MESH
    ////////////////////////////////////////////////////////////////////////

    VertexElementPVec                   pElements;
    VertexConnectionPVec                pConnections;

    uint                              * pVertexPerm;

    ////////////////////////////////////////////////////////////////////////
    // COPIED FROM TETMESH
    ////////////////////////////////////////////////////////////////////////

    uint                                pNTri;
    uint                                pNTet;
    uint                              * pTetrahedrons;
    uint                              * pTriangles;

    TetStubSet                          pTetLUT;

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

STEPS_EXTERN
std::ostream & operator<< (std::ostream & os, steps::solver::efield::TetStub const &);

////////////////////////////////////////////////////////////////////////////////

#endif
// STEPS_SOLVER_EFIELD_TETMESH_HPP

// END
