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

#include "fwd.hpp"

namespace steps::model {

/// Diffusion rule on a vesicle surface a volume system.
///
///\warning Methods start with an underscore are not exposed to Python.

class VesSDiff {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the diffusion rule.
    /// \param vessurfsys Vesicle surface system which the diffusion rule belongs
    /// to. \param lig Reference to the species which the diffusion applies to.
    /// \param dcst Diffusion coefficient of the diffusion rule.
    VesSDiff(std::string const& id, VesSurfsys& vessurfsys, Spec& lig, double dcst = 0.0);

    VesSDiff(const VesSDiff&) = delete;
    VesSDiff& operator=(const VesSDiff&) = delete;

    /// Destructor
    ~VesSDiff();

    ////////////////////////////////////////////////////////////////////////
    // DIFFUSION RULE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the diffusion rule ID.
    ///
    /// \return ID of the diffusion rule.
    inline const std::string& getID() const noexcept {
        return pID;
    }

    /// Set the ID of the diffusion rule.
    ///
    /// \param id ID of the diffusion rule.
    void setID(std::string const& id);

    /// Return a reference to the parent surface system.
    ///
    /// \return Reference to the parent surface system.
    inline VesSurfsys& getVesSurfsys() const noexcept {
        return pVesSurfsys;
    }

    /// Return a reference to the parent model.
    ///
    /// \return Reference to the parent model.
    inline Model& getModel() const noexcept {
        return pModel;
    }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON):
    ////////////////////////////////////////////////////////////////////////

    /// Return a reference to the species to which this diffusion rule applies
    ///
    /// \return Reference of the species
    inline Spec& getLig() const noexcept {
        return *pLig;
    }

    /// Set the species which this difusion rule applies to.
    ///
    /// \param lig Reference to the species
    void setLig(Spec& lig);

    /// Get the rate constant of the diffusion rule.
    ///
    /// \return Rate constant of the diffusion rule.
    inline double getDcst() const noexcept {
        return pDcst;
    }

    /// Set the rate constant of the diffusion rule.
    ///
    /// \param Rate constant of the diffusion rule.
    void setDcst(double dcst);

    ///  Return a list of all species in this diffusion rule.
    ///
    /// \return List of pointers of species.
    /// \warning Currently will return only one species.
    std::vector<Spec*> getAllSpecs() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    Model& pModel;

    VesSurfsys& pVesSurfsys;

    Spec* pLig;
    double pDcst;
};

inline bool operator<(const VesSDiff& lhs, const VesSDiff& rhs) {
    return lhs.getID() < rhs.getID();
}

}  // namespace steps::model
