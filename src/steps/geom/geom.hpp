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
#include <string>
#include <vector>

#include "fwd.hpp"
#include "solver/fwd.hpp"

namespace steps::wm {

/////////////////////////////////////////////////////////////////////////////
/// Geometry container for compartments and patches.
///
/// \warning Methods start with an underscore are not exposed to Python.
////////////////////////////////////////////////////////////////////////////
class Geom {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Geom() = default;
    Geom(const Geom&) = delete;
    Geom& operator=(const Geom&) = delete;

    /// Destructor
    virtual ~Geom();

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON) : COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return a compartment with name id.
    ///
    /// \param id ID of the compartment object.
    /// \return Reference to the compartment.
    wm::Comp& getComp(std::string const& id) const;

    /// Delete a compartment with name id.
    ///
    /// \param id ID of the compartment.
    void delComp(std::string const& id) const;

    /// Return all compartments in the geometry container.
    ///
    /// \return List of pointers to the compartment objects.
    std::vector<wm::Comp*> getAllComps() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return a patch with name id.
    ///
    /// \param id ID of the patch.
    /// \return Reference to the patch.
    wm::Patch& getPatch(std::string const& id) const;

    /// Delete a patch with name id.
    ///
    /// \param id ID of the patch.
    void delPatch(std::string const& id) const;

    /// Return all patches in the geometry container.
    ///
    /// \return List of pointers to the patch objects.
    std::vector<wm::Patch*> getAllPatches() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the compartments in the geometry container.
    ///
    /// \return Number of compartments.
    inline std::size_t _countComps() const noexcept {
        return pComps.size();
    }

    /// Return a compartment with index gidx.
    ///
    /// \param gidx Index of the compartment.
    /// \return Reference to the compartment.
    wm::Comp& _getComp(solver::comp_global_id gidx) const;

    /// Count the patches in the geometry container.
    ///
    /// \return Number of patches.
    inline std::size_t _countPatches() const noexcept {
        return pPatches.size();
    }

    /// Return a patch with index gidx.
    ///
    /// \param gidx Index of the patch.
    /// \return Reference to the patch.
    wm::Patch& _getPatch(solver::patch_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::WM OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a compartment id is occupied.
    ///
    /// \param id ID of the compartment.
    void _checkCompID(std::string const& id) const;

    /// Change the id of a compartment.
    ///
    /// \param o Old id of the compartment.
    /// \param n New id of the compartment.
    void _handleCompIDChange(std::string const& o, std::string const& n);

    /// Add a compartment.
    ///
    /// \param comp Reference to the compartment.
    void _handleCompAdd(wm::Comp& comp);

    /// Delete a compartment.
    ///
    /// \param comp Reference to the compartment.
    void _handleCompDel(wm::Comp& comp);

    /// Check if a patch id is occupied.
    ///
    /// \param id ID of the patch.
    void _checkPatchID(std::string const& id) const;

    /// Change the id of a patch.
    ///
    /// \param o Old id of the patch.
    /// \param n New id of the patch.
    void _handlePatchIDChange(std::string const& o, std::string const& n);

    /// Add a patch.
    ///
    /// \param patch Reference to the patch.
    void _handlePatchAdd(wm::Patch& patch);

    /// Delete a patch.
    ///
    /// \param patch Reference to the patch.
    void _handlePatchDel(wm::Patch& patch);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    std::map<std::string, wm::Comp*> pComps;
    std::map<std::string, wm::Patch*> pPatches;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::wm
