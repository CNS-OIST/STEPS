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
 
#ifndef STEPS_WM_GEOM_HPP
#define STEPS_WM_GEOM_HPP 1


// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace wm {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;
class Patch;

/////////////////////////////////////////////////////////////////////////////
/// Geometry container for compartments and patches.
///
/// \warning Methods start with an underscore are not exposed to Python.
////////////////////////////////////////////////////////////////////////////
class Geom
{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    Geom() = default;

    /// Destructor
    virtual ~Geom();

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON) : COMPARTMENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return a compartment with name id.
    ///
    /// \param id ID of the compartment object.
    /// \return Pointer to the compartment.
    steps::wm::Comp * getComp(std::string const & id) const;

    /// Delete a compartment with name id.
    ///
    /// \param id ID of the compartment.
    void delComp(std::string const & id);

    /// Return all compartments in the geometry container.
    ///
    /// \return List of pointers to the compartment objects.
    std::vector<steps::wm::Comp *> getAllComps() const;

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): PATCHES
    ////////////////////////////////////////////////////////////////////////

    /// Return a patch with name id.
    ///
    /// \param id ID of the patch.
    /// \return Pointer to the patch.
    steps::wm::Patch * getPatch(std::string const & id) const;

    /// Delete a patch with name id.
    ///
    /// \param id ID of the patch.
    void delPatch(std::string const & id);

    /// Return all patches in the geometry container.
    ///
    /// \return List of pointers to the patch objects.
    std::vector<steps::wm::Patch *> getAllPatches() const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Count the compartments in the geometry container.
    ///
    /// \return Number of compartments.
    inline uint _countComps() const
    { return pComps.size(); }

    /// Return a compartment with index gidx.
    ///
    /// \param gidx Index of the compartment.
    /// \return Pointer to the compartment.
    steps::wm::Comp * _getComp(uint gidx) const;

    /// Count the patches in the geometry container.
    ///
    /// \return Number of patches.
    inline uint _countPatches() const
    { return pPatches.size(); }

    /// Return a patch with index gidx.
    ///
    /// \param gidx Index of the patch.
    /// \return Pointer to the patch.
    steps::wm::Patch * _getPatch(uint gidx) const;


    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::WM OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Check if a compartment id is occupied.
    ///
    /// \param id ID of the compartment.
    void _checkCompID(std::string const & id) const;

    /// Change the id of a compartment.
    ///
    /// \param o Old id of the compartment.
    /// \param n New id of the compartment.
    void _handleCompIDChange(std::string const & o, std::string const & n);

    /// Add a compartment.
    ///
    /// \param comp Pointer to the compartment.
    void _handleCompAdd(steps::wm::Comp * comp);

    /// Delete a compartment.
    ///
    /// \param comp Pointer to the compartment.
    void _handleCompDel(steps::wm::Comp * comp);

    /// Check if a patch id is occupied.
    ///
    /// \param id ID of the patch.
    void _checkPatchID(std::string const & id) const;

    /// Change the id of a patch.
    ///
    /// \param o Old id of the patch.
    /// \param n New id of the patch.
    void _handlePatchIDChange(std::string const & o, std::string const & n);

    /// Add a patch.
    ///
    /// \param patch Pointer to the patch.
    void _handlePatchAdd(steps::wm::Patch * patch);

    /// Delete a patch.
    ///
    /// \param patch Pointer to the patch.
    void _handlePatchDel(steps::wm::Patch *patch);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::map<std::string, steps::wm::Comp *>       pComps;
    std::map<std::string, steps::wm::Patch *>      pPatches;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_WM_GEOM_HPP

// END
