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

#include "fwd.hpp"

namespace steps::model {

/// LinkSpecies reactant.
/// Component that represents a reactant that can be referred to from
/// volume and surface systems.
///
/// \warning Methods start with an underscore are not exposed to Python.

class LinkSpec {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the link species.
    /// \param model Reference to the parent model.
    /// \param dcst Default surface diffusion coefficient of the link species
    LinkSpec(std::string const& id, Model& model, double dcst = 0.0);

    LinkSpec(const LinkSpec&) = delete;
    LinkSpec& operator=(const LinkSpec&) = delete;

    /// Destructor
    virtual ~LinkSpec();

    ////////////////////////////////////////////////////////////////////////
    // LINKSPECIES PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the link species ID.
    ///
    /// \return ID of the species.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the link species ID.
    ///
    /// \param id ID of the species.
    void setID(std::string const& id);

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    inline double getDcst() const noexcept {
        return pDcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;

    double pDcst;
};

inline bool operator<(const LinkSpec& lhs, const LinkSpec& rhs) {
    return lhs.getID() < rhs.getID();
}

}  // namespace steps::model
