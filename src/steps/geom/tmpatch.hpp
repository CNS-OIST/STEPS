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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#pragma once

#include "comp.hpp"
#include "patch.hpp"
#include "tetmesh.hpp"
#include "tmcomp.hpp"

#include "math/bbox.hpp"
#include "model/surfsys.hpp"
#include "util/common.h"

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
         std::vector<index_t> const & tris, steps::wm::Comp* wmicomp,
            steps::wm::Comp* wmocomp = nullptr);

    /// Destructor.
    ~TmPatch() {}

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return whether triangles (specified by index) are inside this patch.
    ///
    /// \param tris List of indices of triangles.
    /// \return Results of whether the triangles are inside the patch.
    std::vector<bool> isTriInside(const std::vector<index_t> &tris) const;

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    std::vector<index_t> getAllTriIndices() const
    { return strong_type_to_value_type(pTri_indices); }

    /// Get the minimal coordinate of the rectangular bounding box or a plane.
    ///
    /// \return Minimal coordinate of the rectangular bounding box or a plane.
    std::vector<double> getBoundMin() const;

    /// Get the maximal coordinate of the rectangular bounding box or a plane.
    ///
    /// \return Maximal coordinate of the rectangular bounding box or a plane.
    std::vector<double> getBoundMax() const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all triangles (by index) in the patch.
    ///
    /// \return List of indices of triangles.
    inline std::vector<triangle_id_t> const & _getAllTriIndices() const
    { return pTri_indices; }

    ////////////////////////////////////////////////////////////////////////

private:

    Tetmesh                           * pTetmesh;
    std::vector<triangle_id_t>          pTri_indices;
    std::size_t                         pTrisN;

    steps::math::bounding_box           pBBox;
};

} // namespace tetmesh
} // namespace steps
