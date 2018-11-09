////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2014 Okinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006 University of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPS is free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPS is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program. If not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_TETMESH_SDIFFBOUNDARY_HPP
#define STEPS_TETMESH_SDIFFBOUNDARY_HPP 1

// STL headers
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/geom/sdiffboundary.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/model/surfsys.hpp"
#include "steps/geom/tmpatch.hpp"
#include "steps/geom/geom.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;
class SDiffBoundary;

// Auxiliary declarations.
typedef SDiffBoundary *                          SDiffBoundaryP;
typedef std::map<std::string, SDiffBoundaryP>    SDiffBoundaryPMap;
typedef SDiffBoundaryPMap::iterator              SDiffBoundaryPMapI;
typedef SDiffBoundaryPMap::const_iterator        SDiffBoundaryPMapCI;

typedef std::vector<SDiffBoundaryP>              SDiffBoundaryPVec;
typedef SDiffBoundaryPVec::iterator              SDiffBoundaryPVecI;
typedef SDiffBoundaryPVec::const_iterator        SDiffBoundaryPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of surface diffusion boundary bars of a Tetmesh.
///
/// Tetmesh object is responsible for maintaining lifetime of associated
/// SDiffBoundary objects (Python proxy class must set thisown to zero.)
///
/// \warning Methods start with an underscore are not exposed to Python.
class SDiffBoundary
{

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
    SDiffBoundary(std::string id, Tetmesh * container,
            std::vector<uint> const & bars, std::vector<steps::tetmesh::TmPatch *> const & patches);

    /// Destructor.
    virtual ~SDiffBoundary() {}

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return the patch id.
    ///
    /// \return ID of the surface diffusion boundary.
    std::string const & getID() const
    { return pID; }

    /// Set or change the surface diffusion boundary id.
    ///
    /// \param id ID of the surface diffusion boundary.
    void setID(std::string const & id);

    /// Return a pointer to the geometry container object.
    ///
    /// \return Pointer to the parent geometry container.
    steps::tetmesh::Tetmesh * getContainer() const
    { return pTetmesh; }

    /// Return whether bars (specified by index) are inside this surface diffusion boundary.
    ///
    /// \param bar List of indices of bars.
    /// \return Results of whether the bars are inside the surface diffusion boundary.
    std::vector<bool> isBarInside(const std::vector<uint> &bars) const;

    /// Return all bars (by index) in the surface diffusion boundary.
    ///
    /// \return List of indices of bars.
    inline std::vector<uint> const & getAllBarIndices() const
    { return pBar_indices; }

    /// Return the patches this surface diffusion boundary connects
    ///
    /// \return List of the two patches.
    inline std::vector<steps::wm::Patch *> getPatches() const
    { return std::vector<steps::wm::Patch *>{pIPatch, pOPatch}; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all bars (by index) in the surface diffusion boundary.
    ///
    /// \return List of indices of bars.
    inline std::vector<uint> const & _getAllBarIndices() const
    { return pBar_indices; }

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;

    steps::wm::Patch *                  pIPatch;
    steps::wm::Patch *                  pOPatch;

    steps::tetmesh::Tetmesh           * pTetmesh;
    std::vector<uint>                   pBar_indices;
    uint                                pBarsN;

////////////////////////////////////////////////////////////////////////////////

};

}
}

#endif
// STEPS_TETMESH_SDIFFBOUNDARY_HPP

// END
