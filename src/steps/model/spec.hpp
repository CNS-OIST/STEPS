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

/// Species reactant.
/// Component that represents a reactant that can be referred to from
/// volume and surface systems.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Spec {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the species.
    /// \param model Reference to the parent model.
    Spec(std::string const& id, Model& model, int valence = 0);

    Spec(const Spec&) = delete;
    Spec& operator=(const Spec&) = delete;

    /// Destructor
    virtual ~Spec();

    ////////////////////////////////////////////////////////////////////////
    // SPECIES PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the species ID.
    ///
    /// \return ID of the species.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the species ID.
    ///
    /// \param id ID of the species.
    virtual void setID(std::string const& id);

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    /// Set the valence of the species.
    ///
    /// \param valence Valence of the species.
    void setValence(int valence);

    /// Return the valence of the species.
    ///
    /// \return Valence of the species.
    inline int getValence() const noexcept {
        return pValence;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    virtual void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;
    int pValence;
};

inline bool operator<(const Spec& lhs, const Spec& rhs) {
    return lhs.getID() < rhs.getID();
}

}  // namespace steps::model
