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

#include <map>
#include <set>
#include <string>
#include <vector>

#include "fwd.hpp"

#include "model/fwd.hpp"
#include "util/collections.hpp"

namespace steps::wm {

/// Base class for compartment objects.
///
///    It provides basic functionality and data that is shared by all classes
///    derived from Comp:
///        - Getting and setting a valid compartment ID string, and handling
///          the interaction with the container object.
///        - Getting (and at least in this base class also setting) the total
///          volume of the compartment.
///        - The volume systems implemented in the compartment.
///        - References to Patch objects.
///
///    This base class can be used directly with well-mixed solvers.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Comp {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the compartment.
    /// \param container Reference to the parent geometry container.
    /// \param vol Volume of the compartment.
    Comp(std::string id, wm::Geom& container, double vol = 0.0);

    Comp(const Comp&) = delete;
    Comp& operator=(const Comp&) = delete;

    /// Destructor
    virtual ~Comp();

    ////////////////////////////////////////////////////////////////////////
    // COMPARTMENT PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the compartment id.
    ///
    /// \return ID of the compartment.
    inline std::string const& getID() const noexcept {
        return pID;
    }

    /// Set or change the compartment id.
    ///
    /// \param id ID of the compartment.
    void setID(std::string const& id);

    /// Return a reference to the geometry container object.
    ///
    /// \return Reference to the parent geometry container.
    inline wm::Geom& getContainer() const noexcept {
        return pContainer;
    }

    /// Return the volume of the compartment.
    ///
    /// \return Volume of the compartment.
    inline double getVol() const noexcept {
        return pVol;
    }

    /// Set the volume of the compartment.
    ///
    /// \param vol Volume of the compartment.
    virtual void setVol(double vol);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): VOLUME SYSTEM
    ////////////////////////////////////////////////////////////////////////

    /// Add a volume system with name id.
    ///
    /// \param id ID of the volume system.
    void addVolsys(std::string const& id);

    /// Return a list of volume systems implemented by the compartment.
    ///
    /// \return List of ids of volume systems.
    inline const auto& getVolsys() const noexcept {
        return pVolsys;
    }

    /// Delete a volume system with name id.
    ///
    /// \param id ID of the volume system.
    void delVolsys(std::string const& id);

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): MODEL LINKING
    ////////////////////////////////////////////////////////////////////////
    /// Return all spec in the compartment giving a model.
    util::flat_set<model::Spec*> getAllSpecs(const model::Model& model) const;

    /// Return all Reac in the compartment giving a model.
    util::flat_set<model::Reac*> getAllReacs(const model::Model& model) const;

    /// Return all Diff in the compartment giving a model.
    util::flat_set<model::Diff*> getAllDiffs(const model::Model& model) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS (EXPOSED TO PYTHON): PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return the set of inner patches.
    ///
    /// \return List of pointers to the inner patches.
    /// \warning If a compartment is set as the outer compartment of a patch,
    ///          this patch is a inner patch of the compartment.
    inline const auto& getIPatches() const noexcept {
        return pIPatches;
    }

    /// Return the set of outer patches.
    ///
    /// \return List of pointers to the outer patches.
    /// \warning If a compartment is set as the inner compartment of a patch,
    ///          this patch is a outer patch of the compartment.
    inline const auto& getOPatches() const noexcept {
        return pOPatches;
    }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Add a inner patch.
    ///
    /// \param patch Reference to the inner patch.
    void _addIPatch(wm::Patch& patch);

    /// Delete an inner patch.
    ///
    /// \param patch Reference to the inner patch.
    void _delIPatch(wm::Patch& patch);

    /// Add an outer patch.
    ///
    /// \param patch Reference to the outer patch.
    void _addOPatch(wm::Patch& patch);

    /// Delete an outer patch.
    ///
    /// \param patch Reference to the outer patch.
    void _delOPatch(wm::Patch& patch);

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

    /// Count the inner patches.
    ///
    /// \return Number of the inner patches.
    inline uint _countIPatches() const noexcept {
        return pIPatches.size();
    }

    /// Get a inner patch with index lidx.
    ///
    /// \param lidx Index of the inner patch.
    /// \return Pointer to the inner patch.
    wm::Patch* _getIPatch(uint lidx) const;

    /// Count the outer patches.
    ///
    /// \return Number of the outer patches.
    inline uint _countOPatches() const noexcept {
        return pOPatches.size();
    }

    /// Get a outer patch with index lidx.
    ///
    /// \param lidx Index of the outer patch.
    /// \return Pointer to the outer patch.
    wm::Patch* _getOPatch(uint lidx) const;

    ////////////////////////////////////////////////////////////////////////

  protected:
    ////////////////////////////////////////////////////////////////////////

    double pVol;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::string pID;
    wm::Geom& pContainer;
    util::flat_set<std::string> pVolsys;
    util::flat_set<wm::Patch*> pIPatches;
    util::flat_set<wm::Patch*> pOPatches;
};

}  // namespace steps::wm
