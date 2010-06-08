////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_WM_COMP_HPP
#define STEPS_WM_COMP_HPP 1


// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <set>

// STEPS headers.
#include "../common.h"
#include "../model/volsys.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wm)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Geom;
class Patch;
class Comp;

// Auxiliary declarations.
typedef Comp *                          CompP;
typedef std::map<std::string, CompP>    CompPMap;
typedef CompPMap::iterator              CompPMapI;
typedef CompPMap::const_iterator        CompPMapCI;

typedef std::vector<CompP>              CompPVec;
typedef CompPVec::iterator              CompPVecI;
typedef CompPVec::const_iterator        CompPVecCI;

////////////////////////////////////////////////////////////////////////////////

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

class Comp
{
public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the compartment.
    /// \param container Pointer to the parent geometry container.
    /// \param volsys Pointer to the implemented volume system.
    /// \param vol Volume of the compartment.
	Comp(std::string const & id, steps::wm::Geom * container, double vol = 0.0);

    /// Destructor
	virtual ~Comp(void);

	////////////////////////////////////////////////////////////////////////
	// COMPARTMENT PROPERTIES
	////////////////////////////////////////////////////////////////////////

	/// Return the compartment id.
    ///
    /// \return ID of the compartment.
	std::string getID(void) const
	{ return pID; }

	/// Set or change the compartment id.
    ///
    /// \param id ID of the compartment.
	void setID(std::string const & id);

	/// Return a pointer to the geometry container object.
    ///
    /// \return Pointer to the parent geometry container.
	steps::wm::Geom * getContainer(void) const
	{ return pContainer; }

    /// Return the volume of the compartment.
    ///
    /// \return Volume of the compartment.
	double getVol(void) const
	{ return pVol; }

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
	void addVolsys(std::string const & id);

    /// Return a list of volume systems implemented by the compartment.
    ///
    /// \return List of ids of volume systems.
	std::set<std::string> getVolsys(void) const
	{ return pVolsys; }

    /// Delete a volume system with name id.
    ///
    /// \param id ID of the volume system.
	void delVolsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): PATCHES
	////////////////////////////////////////////////////////////////////////

	/// Return a copy of the set of inner patches.
    ///
    /// \return List of pointers to the inner patches.
    /// \warning If a compartment is set as the outer compartment of a patch,
    ///          this patch is a inner patch of the compartment.
	std::set<steps::wm::Patch *> getIPatches(void) const
	{ return pIPatches; }

	/// Return a copy of the set of outer patches.
    ///
    /// \return List of pointers to the outer patches.
    /// \warning If a compartment is set as the inner compartment of a patch,
    ///          this patch is a outer patch of the compartment.
	std::set<steps::wm::Patch *> getOPatches(void) const
	{ return pOPatches; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: PATCHES
	////////////////////////////////////////////////////////////////////////

    /// Add a inner patch.
    ///
    /// \param patch Pointer to the inner patch.
	void _addIPatch(steps::wm::Patch * patch);

    /// Delete an inner patch.
    ///
    /// \param patch Pointer to the inner patch.
	void _delIPatch(steps::wm::Patch * patch);

    /// Add an outer patch.
    ///
    /// \param patch Pointer to the outer patch.
	void _addOPatch(steps::wm::Patch * patch);

    /// Delete an outer patch.
    ///
    /// \param patch Pointer to the outer patch.
	void _delOPatch(steps::wm::Patch * patch);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
	/// Called if Python object deleted, or from del method in parent object.
	/// Will only be called once
	void _handleSelfDelete(void);


	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

    /// Count the inner patches.
    ///
    /// \return Number of the inner patches.
	inline uint _countIPatches(void) const
	{ return pIPatches.size(); }

    /// Get a inner patch with index lidx.
    ///
    /// \param lidx Index of the inner patch.
    /// \return Pointer to the inner patch.
	steps::wm::Patch * _getIPatch(uint lidx) const;

    /// Count the outer patches.
    ///
    /// \return Number of the outer patches.
	inline uint _countOPatches(void) const
	{ return pIPatches.size(); }

    /// Get a outer patch with index lidx.
    ///
    /// \param lidx Index of the outer patch.
    /// \return Pointer to the outer patch.
	steps::wm::Patch * _getOPatch(uint lidx) const;

	////////////////////////////////////////////////////////////////////////

protected:

	////////////////////////////////////////////////////////////////////////

	double                              pVol;

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;
	steps::wm::Geom                   * pContainer;
	std::set<std::string>               pVolsys;
	std::set<steps::wm::Patch *>        pIPatches;
	std::set<steps::wm::Patch *>        pOPatches;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wm)
END_NAMESPACE(steps)

#endif
// STEPS_WM_COMP_HPP
