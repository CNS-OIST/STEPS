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

#include "comp.hpp"
#include "math/bbox.hpp"
#include "util/common.hpp"
#include "util/vocabulary.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::tetmesh {

////////////////////////////////////////////////////////////////////////////////

// Forward & auxiliary declarations.
class Tetmesh;

////////////////////////////////////////////////////////////////////////////////

/// Provides annotation for a group of tetrahedron in a Tetmesh.
///
///
/// \warning Methods start with an underscore are not exposed to Python.
class TmComp: public wm::Comp {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor.
    ///
    /// \param id ID of the TmComp.
    /// \param container Temesh container for the tetrahedrons.
    /// \param tets A sequence of tetrahedron (by index) as a vector
    ///             of unsigned integers which is represented as a
    ///             sequence of positive integer values) in Python.
    /// \param volsys Pointer to the volume system associated.
    ///
    TmComp(std::string const& id, Tetmesh* container, std::vector<index_t> const& tets);

    ~TmComp();
    ////////////////////////////////////////////////////////////////////////
    // BASE CLASS METHODS
    ////////////////////////////////////////////////////////////////////////

    void setVol(double vol) override;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a list of all tetrahedron by indices.
    ///
    /// \return List of indices of the tetrahedrons.
    inline const std::vector<index_t> getAllTetIndices() const noexcept {
        return strong_type_to_value_type(pTet_indices);
    }

    /// Return the number of tetrahedrons in this TmComp
    ///
    /// \return the number of tetrahedrons in this TmCOmp
    inline uint countTets() const noexcept {
        return pTetsN;
    }

    // Return whether tetrahedrons (specified by index) are inside this
    // compartment.
    ///
    /// \param tet List of indices of tetrahedrons.
    /// \return List of results of the tetrahedrons are inside the compartment.
    std::vector<bool> isTetInside(const std::vector<index_t>& tets) const;

    /// Get the minimal coordinate of the rectangular bounding box.
    ///
    /// \return Minimal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMin() const;

    /// Get the maximal coordinate of the rectangular bounding box.
    ///
    /// \return Maximal coordinate of the rectangular bounding box.
    std::vector<double> getBoundMax() const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO C++)
    ////////////////////////////////////////////////////////////////////////

    /// Return all tetrahedrons (by index) in the compartment.
    ///
    /// \return List of indices of tetrahedrons.
    inline std::vector<tetrahedron_global_id> const& _getAllTetIndices() const noexcept {
        return pTet_indices;
    }

  private:
    Tetmesh* pTetmesh;
    std::vector<tetrahedron_global_id> pTet_indices;
    std::size_t pTetsN{0};
    math::bounding_box pBBox;
};

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::tetmesh
