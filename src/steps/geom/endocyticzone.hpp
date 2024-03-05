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

#include "fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps::tetmesh {

/// Endocytic zone composed of a group of surface triangles of a TmPatch.
///
/// \warning Methods start with an underscore are not exposed to Python.
class EndocyticZone {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the endocytic zone
    /// \param patch patch in which the zone is defined
    /// \param tris A sequence of triangles (by index)
    EndocyticZone(std::string const& id, TmPatch& patch, std::vector<index_t> const& tris);
    EndocyticZone(const EndocyticZone&) = delete;
    EndocyticZone& operator=(const EndocyticZone&) = delete;

    /// Destructor.
    ~EndocyticZone();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return the endocytic zone id.
    ///
    /// \return ID of the endocytic zone.
    inline std::string const& getID() const noexcept {
        return pID;
    }

    /// Return a reference to the patch container object.
    ///
    /// \return The parent patch
    inline TmPatch& getPatch() const noexcept {
        return pPatch;
    }

    /// Return all triangles (by index) in the endocytic zone.
    ///
    /// \return List of indices of triangles.
    inline std::vector<triangle_global_id> const& getAllTriIndices() const noexcept {
        return pTri_indices;
    }


    ////////////////////////////////////////////////////////////////////////

  private:
    std::string pID;
    TmPatch& pPatch;
    std::vector<triangle_global_id> pTri_indices;
};

}  // namespace steps::tetmesh
