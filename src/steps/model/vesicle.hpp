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
#include "util/collections.hpp"

namespace steps::model {

/// Vesicles .
/// Component that represents a vesicle that can be referred to from
/// volume system.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Vesicle {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the vesicles.
    /// \param model Reference to the parent model.
    /// \param diameter Diameter of the vesicle (m)
    /// \param dcst Default diffusion coefficient of the vesicle (m2/s)

    Vesicle(std::string const& id, Model& model, double diameter, double dcst = 0.0);

    Vesicle(const Vesicle&) = delete;
    Vesicle& operator=(const Vesicle&) = delete;

    /// Destructor
    virtual ~Vesicle();

    ////////////////////////////////////////////////////////////////////////
    // VESICLEIES PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the vesicle's ID.
    ///
    /// \return ID of the vesicles.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set or change the vesicles ID.
    ///
    /// \param id ID of the vesicles.
    void setID(std::string const& id);

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    /// Return the vesicle's diameter.
    ///
    /// \return Diameter of the vesicle.
    inline double getDiameter() const noexcept {
        return pDiameter;
    }

    /// Get the rate constant of the vesicle diffusion.
    ///
    /// \return Rate constant of the vesicle diffusion .
    inline double getDcst() const noexcept {
        return pDcst;
    }

    /// Set the rate constant of the vesicle diffusion.
    ///
    /// \param Rate constant of the vesicle diffusion.
    void setDcst(double dcst);

    /// Associate a vesicle surface system with name id.
    ///
    /// \param id ID of the vesicle surface system.
    void addVesSurfsys(std::string const& id);

    /// Get a vesicle surface system.
    ///
    /// \return List of the vesicle surface systems associated to the vesicle.
    inline const auto& getVesSurfsys() const noexcept {
        return pVesSurfsys;
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
    double pDiameter;

    double pDcst;

    util::flat_set<std::string> pVesSurfsys;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::model
