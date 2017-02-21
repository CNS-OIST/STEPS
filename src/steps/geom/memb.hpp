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

#ifndef STEPS_TETMESH_MEMB_HPP
#define STEPS_TETMESH_MEMB_HPP 1

// STEPS headers.
#include "steps/common.h"
#include "steps/geom/patch.hpp"
//#include "tetmesh.hpp"

// STL headers
#include <vector>
#include <string>
//#include <iostream>

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Tetmesh;
class Memb;
class TmPatch;

// Auxiliary declarations.
typedef Memb *                          MembP;
typedef std::map<std::string, MembP>    MembPMap;
typedef MembPMap::iterator              MembPMapI;
typedef MembPMap::const_iterator        MembPMapCI;

typedef std::vector<MembP>              MembPVec;
typedef MembPVec::iterator              MembPVecI;
typedef MembPVec::const_iterator        MembPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a membrane in a Tetmesh. If membrane potential
/// is included in simulation, the potential will be calculated for this
/// closed group of triangles in the tetrahedral mesh.
/// Currently the only check this object performs is a closed-surface test,
/// all checks on Patches on this membrane will be made during Tetmesh
/// construction
///
/// Tetmesh object is responsible for maintaining lifetime of associated
/// Memb objects (Python proxy class must set thisown to zero.)
///
/// \warning Methods start with an underscore are not exposed to Python.

class Memb
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the TmPatch.
    /// \param container Pointer to the Tetmesh container.
    /// \param patches A sequence of TmPatches as a vector
    ///             of pointers which is represented as
    ///             a sequence in Python.
    ///
    Memb(std::string const & id, Tetmesh * container,
            std::vector<TmPatch *> const & patches,
            bool verify = false, uint opt_method = 1, double search_percent = 100.0, std::string const & opt_file_name = std::string());

    // Destructor
    ~Memb() {}

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the tetmesh container object.
    ///
    /// \return Pointer to the parent tetmesh container.
    steps::tetmesh::Tetmesh * getContainer(void) const
    { return pTetmesh; }

    /// Return the membrane id.
    ///
    /// \return ID of the membrane.
    std::string getID(void) const
    { return pID; }

    /// Return whether triangles (specified by index) are inside this patch.
    ///
    /// \param tri List of indices of triangles.
    /// \return Results of whether the triangles are inside the patch.
    std::vector<bool> isTriInside(const std::vector<uint> &tri) const;

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<uint> getAllTriIndices(void) const
    { return pTri_indices; }

    /// Return the number of triangles in this Memb.
    ///
    /// \return the number of triangles in this Memb
    inline uint countTris(void) const
    { return pTrisN; }

    /// Return all tetrahedrons in the conduction volume.
    ///
    /// \return List of indices of tetrahedrons in conduction volume.
    inline std::vector<uint> getAllVolTetIndices(void) const
    { return pTet_indices; }

    /// Return the number of tetrahedrons in the conduction volume.
    ///
    /// \return the number of tetrahedrons in the conduction volume
    inline uint countVolTets(void) const
    { return pTetsN; }

    /// Return all 'virtual triangles' that is the triangles that
    /// make up the closed surface if supplied triangles form an open surface.
    /// These triangles will be voltage-clamped for the time being with
    /// the possibility of coupling to outside whole-cell simulators in
    /// the future.
    ///
    /// \return List of indices of virtual triangles.
    inline std::vector<uint> getAllVirtTriIndices(void) const
    { return pTrivirt_indices; }

    /// Return the number of virtual triangles in this Memb.
    ///
    /// \return the number of virtual triangles in this Memb
    inline uint countVirtTris(void) const
    { return pTriVirtsN; }

    /// Return all vertices in the conduction volume and membrane surface.
    ///
    /// \return List of indices of vertices in membrane surface and conduction volume.
    inline std::vector<uint> getAllVertIndices(void) const
    { return pVert_indices; }

    /// Return the number of vertices in the conduction volume and membrane surface.
    ///
    /// \return the number of vertices in membrane surface and conduction volume.
    inline uint countVerts(void) const
    { return pVertsN; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<uint> const & _getAllTriIndices(void) const
    { return pTri_indices; }

    /// Return all tetrahedrons in the conduction volume.
    ///
    /// \return List of indices of tetrahedrons in conduction volume.
    inline std::vector<uint> const & _getAllVolTetIndices(void) const
    { return pTet_indices; }

    /// Return all 'virtual triangles' that is the triangles that
    /// make up the closed surface if supplied triangles form an open surface.
    /// These triangles will be voltage-clamped for the time being with
    /// the possibility of coupling to outside whole-cell simulators in
    /// the future.
    ///
    /// \return List of indices of virtual triangles.
    inline std::vector<uint> const & _getAllVirtTriIndices(void) const
    { return pTrivirt_indices; }

    /// Return all vertices in the conduction volume and membrane surface.
    ///
    /// \return List of indices of vertices in membrane surface and conduction volume.
    inline std::vector<uint> const & _getAllVertIndices(void) const
    { return pVert_indices; }

    /// Return whether surface is 'open' or not
    ///
    /// \return Bool of open or not.
    inline bool open(void) const
    { return pOpen; }

    /// Return the method the user has specified to optimize the vertex
    /// indexing in the conduction volume.
    inline uint _getOpt_method(void)
    { return pOpt_method; }

    /// Return the percentage of starting nodes tested for breadth-first search (default 100%)
    /// Will be ignored if principal axis search is used
    inline double _getSearch_percent(void)
    { return pSearch_percent; }


    /// Return optimization file, if it exists (otherwise empty string)
    std::string _getOpt_file_name(void) const
    { return pOpt_file_name; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    void verifyMemb();

    std::string                         pID;
    Tetmesh                           * pTetmesh;

    std::vector<uint>                   pTri_indices;

    // The conduction volume tetrahedron indices
    std::vector<uint>                   pTet_indices;

    // The 'virtual triangles', triangles that will have to be voltage-clamped
    std::vector<uint>                   pTrivirt_indices;

    // The vertices
    std::vector<uint>                   pVert_indices;

    uint                                pTrisN;
    uint                                pTetsN;
    uint                                pTriVirtsN;
    uint                                pVertsN;

    bool                                pOpen;

    uint                                pOpt_method;

    std::string                         pOpt_file_name;

    double                              pSearch_percent;

////////////////////////////////////////////////////////////////////////////////

};

}
}

#endif

// STEPS_TETMESH_MEMB

// END
