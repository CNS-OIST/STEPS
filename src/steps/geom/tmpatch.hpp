/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_TETMESH_TMPATCH_HPP
#define STEPS_TETMESH_TMPATCH_HPP 1


// STEPS headers.
#include "steps/common.h"
#include "steps/geom/patch.hpp"
#include "steps/geom/tetmesh.hpp"
#include "steps/geom/tmcomp.hpp"
#include "steps/geom/comp.hpp"
#include "steps/math/bbox.hpp"
#include "steps/model/surfsys.hpp"

// STL headers
#include <vector>

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;
class TmComp;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of surface triangles of a Tetmesh.
///
/// \warning Methods start with an underscore are not exposed to Python.
class TmPatch : public steps::wm::Patch
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the TmPatch.
    /// \param container Pointer to the Tetmesh container.
    /// \param tris A sequence of triangles (by index) as a vector
    ///             of unsigned integers which is represented as
    ///             a sequence of positive integer values) in Python.
    /// \param icomp Pointer to the inner compartment.
    /// \param ocomp Pointer to the outer compartment.
    /// \param surfsys Pointer to the assocaited surface system.
    ///
    /// This is the constructor for the wm (well-mixed) namespace.
    TmPatch(std::string const & id, Tetmesh * container,
         std::vector<uint> const & tris, steps::wm::Comp* icomp,
            steps::wm::Comp* ocomp = 0);

    /// Destructor.
    ~TmPatch() {}

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return whether triangles (specified by index) are inside this patch.
    ///
    /// \param tris List of indices of triangles.
    /// \return Results of whether the triangles are inside the patch.
    std::vector<bool> isTriInside(const std::vector<uint> &tris) const;

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    std::vector<uint> getAllTriIndices(void) const
    { return pTri_indices; }

    /// Get the minimal coordinate of the rectangular bounding box or a plane.
    ///
    /// \return Minimal coordinate of the rectangular bounding box or a plane.
    std::vector<double> getBoundMin(void) const;

    /// Get the maximal coordinate of the rectangular bounding box or a plane.
    ///
    /// \return Maximal coordinate of the rectangular bounding box or a plane.
    std::vector<double> getBoundMax(void) const;
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<uint> const & _getAllTriIndices(void) const
    { return pTri_indices; }

    ////////////////////////////////////////////////////////////////////////

private:

    Tetmesh                           * pTetmesh;
    std::vector<uint>                   pTri_indices;
    uint                                pTrisN;
    
    steps::math::bounding_box           pBBox;
};

}
}

#endif
// STEPS_TETMESH_TMPATCH_HPP

// END
