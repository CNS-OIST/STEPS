/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

#include <string>
#include <vector>

#include "patch.hpp"
#include "tmcomp.hpp"

#include "util/common.hpp"
#include "util/vocabulary.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Auxiliary declarations.
typedef Memb* MembP;
typedef std::map<std::string, MembP> MembPMap;
typedef MembPMap::iterator MembPMapI;
typedef MembPMap::const_iterator MembPMapCI;

typedef std::vector<MembP> MembPVec;
typedef MembPVec::iterator MembPVecI;
typedef MembPVec::const_iterator MembPVecCI;

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

class Memb {
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
    Memb(std::string id,
         Tetmesh* container,
         std::vector<TmPatch*> const& patches,
         std::vector<TmComp*> const& compartments,
         bool verify = false,
         uint opt_method = 1,
         double search_percent = 100.0,
         std::string opt_file_name = {});

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the tetmesh container object.
    ///
    /// \return Pointer to the parent tetmesh container.
    inline tetmesh::Tetmesh* getContainer() const noexcept {
        return pTetmesh;
    }

    /// Return the membrane id.
    ///
    /// \return ID of the membrane.
    inline std::string const& getID() const noexcept {
        return pID;
    }

    /// Return whether triangles (specified by index) are inside this patch.
    ///
    /// \param tri List of indices of triangles.
    /// \return Results of whether the triangles are inside the patch.
    std::vector<bool> isTriInside(const std::vector<index_t>& tri) const;

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<index_t> getAllTriIndices() const noexcept {
        return strong_type_to_value_type(pTri_indices);
    }

    /// Return the number of triangles in this Memb.
    ///
    /// \return the number of triangles in this Memb
    inline uint countTris() const noexcept {
        return pTrisN;
    }

    /// Return all tetrahedrons in the conduction volume.
    ///
    /// \return List of indices of tetrahedrons in conduction volume.
    inline std::vector<index_t> getAllVolTetIndices() const noexcept {
        return strong_type_to_value_type(pTet_indices);
    }

    /// Return the number of tetrahedrons in the conduction volume.
    ///
    /// \return the number of tetrahedrons in the conduction volume
    inline uint countVolTets() const noexcept {
        return pTetsN;
    }

    /// Return all 'virtual triangles' that is the triangles that
    /// make up the closed surface if supplied triangles form an open surface.
    /// These triangles will be voltage-clamped for the time being with
    /// the possibility of coupling to outside whole-cell simulators in
    /// the future.
    ///
    /// \return List of indices of virtual triangles.
    inline std::vector<index_t> getAllVirtTriIndices() const noexcept {
        return strong_type_to_value_type(pTrivirt_indices);
    }

    /// Return the number of virtual triangles in this Memb.
    ///
    /// \return the number of virtual triangles in this Memb
    inline uint countVirtTris() const noexcept {
        return pTriVirtsN;
    }

    /// Return all vertices in the conduction volume and membrane surface.
    ///
    /// \return List of indices of vertices in membrane surface and conduction
    /// volume.
    inline std::vector<index_t> getAllVertIndices() const noexcept {
        return strong_type_to_value_type(pVert_indices);
    }

    /// Return the number of vertices in the conduction volume and membrane
    /// surface.
    ///
    /// \return the number of vertices in membrane surface and conduction volume.
    inline uint countVerts() const noexcept {
        return pVertsN;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<triangle_global_id> const& _getAllTriIndices() const noexcept {
        return pTri_indices;
    }

    /// Return all tetrahedrons in the conduction volume.
    ///
    /// \return List of indices of tetrahedrons in conduction volume.
    inline std::vector<tetrahedron_global_id> const& _getAllVolTetIndices() const noexcept {
        return pTet_indices;
    }

    /// Return all 'virtual triangles' that is the triangles that
    /// make up the closed surface if supplied triangles form an open surface.
    /// These triangles will be voltage-clamped for the time being with
    /// the possibility of coupling to outside whole-cell simulators in
    /// the future.
    ///
    /// \return List of indices of virtual triangles.
    inline std::vector<triangle_global_id> const& _getAllVirtTriIndices() const noexcept {
        return pTrivirt_indices;
    }

    /// Return all vertices in the conduction volume and membrane surface.
    ///
    /// \return List of indices of vertices in membrane surface and conduction
    /// volume.
    inline std::vector<vertex_id_t> const& _getAllVertIndices() const noexcept {
        return pVert_indices;
    }

    /// Return whether surface is 'open' or not
    ///
    /// \return Bool of open or not.
    inline bool open() const noexcept {
        return pOpen;
    }

    /// Return the method the user has specified to optimize the vertex
    /// indexing in the conduction volume.
    inline uint _getOpt_method() const noexcept {
        return pOpt_method;
    }

    /// Return the percentage of starting nodes tested for breadth-first search
    /// (default 100%) Will be ignored if principal axis search is used
    inline double _getSearch_percent() const noexcept {
        return pSearch_percent;
    }

    /// Return optimization file, if it exists (otherwise empty string)
    inline std::string const& _getOpt_file_name() const noexcept {
        return pOpt_file_name;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    void verifyMemb();

    std::string pID;
    Tetmesh* pTetmesh;

    std::vector<triangle_global_id> pTri_indices;

    // The conduction volume tetrahedron indices
    std::vector<tetrahedron_global_id> pTet_indices;

    // The 'virtual triangles', triangles that will have to be voltage-clamped
    std::vector<triangle_global_id> pTrivirt_indices;

    // The vertices
    std::vector<vertex_id_t> pVert_indices;

    uint pTrisN{0};
    uint pTetsN{0};
    uint pTriVirtsN{0};
    uint pVertsN{0};

    bool pOpen{false};

    uint pOpt_method;

    std::string pOpt_file_name;

    double pSearch_percent;

    ////////////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::tetmesh
