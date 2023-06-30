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

#include "geom.hpp"
#include "tetmesh.hpp"
#include "tmpatch.hpp"

#include "model/surfsys.hpp"
#include "util/common.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;
class SDiffBoundary;

// Auxiliary declarations.
typedef SDiffBoundary* SDiffBoundaryP;
typedef std::map<std::string, SDiffBoundaryP> SDiffBoundaryPMap;
typedef SDiffBoundaryPMap::iterator SDiffBoundaryPMapI;
typedef SDiffBoundaryPMap::const_iterator SDiffBoundaryPMapCI;

typedef std::vector<SDiffBoundaryP> SDiffBoundaryPVec;
typedef SDiffBoundaryPVec::iterator SDiffBoundaryPVecI;
typedef SDiffBoundaryPVec::const_iterator SDiffBoundaryPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of surface diffusion boundary bars of a
/// Tetmesh.
///
/// Tetmesh object is responsible for maintaining lifetime of associated
/// SDiffBoundary objects (Python proxy class must set thisown to zero.)
///
/// \warning Methods starting with an underscore are not exposed to Python.
class SDiffBoundary {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the SDiffBoundary.
    /// \param container Pointer to the Tetmesh container.
    /// \param bars A sequence of bars (by index) as a vector
    ///             of unsigned integers which is represented as
    ///             a sequence of positive integer values in Python.
    /// \param patches A sequence of patches (length 2) as a vector
    ///             of steps.tetmesh.TmPatch pointers which is represented as
    ///             a sequence of TmPatch object references in Python.

    ///
    /// This is the constructor for the tetmesh (tetrahedron mesh) namespace.
    SDiffBoundary(std::string id,
                  Tetmesh* container,
                  std::vector<index_t> const& bars,
                  std::vector<tetmesh::TmPatch*> const& patches);

    /// Destructor.
    virtual ~SDiffBoundary() {}

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return the patch id.
    ///
    /// \return ID of the surface diffusion boundary.
    inline std::string const& getID() const noexcept {
        return pID;
    }

    /// Set or change the surface diffusion boundary id.
    ///
    /// \param id ID of the surface diffusion boundary.
    void setID(std::string const& id);

    /// Return a pointer to the geometry container object.
    ///
    /// \return Pointer to the parent geometry container.
    inline tetmesh::Tetmesh* getContainer() const noexcept {
        return pTetmesh;
    }

    /// Return whether bars (specified by index) are inside this surface diffusion
    /// boundary.
    ///
    /// \param bar List of indices of bars.
    /// \return Results of whether the bars are inside the surface diffusion
    /// boundary.
    std::vector<bool> isBarInside(const std::vector<index_t>& bars) const;

    /// Return all bars (by index) in the surface diffusion boundary.
    ///
    /// \return List of indices of bars.
    inline std::vector<index_t> const& getAllBarIndices() const noexcept {
        return pBar_indices;
    }

    /// Return the patches this surface diffusion boundary connects
    ///
    /// \return List of the two patches.
    inline std::vector<wm::Patch*> getPatches() const noexcept {
        return {pIPatch, pOPatch};
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all bars (by index) in the surface diffusion boundary.
    ///
    /// \return List of indices of bars.
    inline std::vector<index_t> const& _getAllBarIndices() const noexcept {
        return pBar_indices;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;

    wm::Patch* pIPatch{nullptr};
    wm::Patch* pOPatch{nullptr};

    tetmesh::Tetmesh* pTetmesh;
    std::vector<index_t> pBar_indices;
};

}  // namespace steps::tetmesh
