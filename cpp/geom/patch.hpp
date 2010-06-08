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

#ifndef STEPS_WM_PATCH_HPP
#define STEPS_WM_PATCH_HPP 1


// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <set>

// STEPS headers.
#include "../common.h"
#include "../model/surfsys.hpp"
#include "geom.hpp"
#include "comp.hpp"

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wm)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.

class Comp;
class Geom;
class Patch;

// Auxiliary declarations.
typedef Patch *                         PatchP;
typedef std::map<std::string, PatchP>   PatchPMap;
typedef PatchPMap::iterator             PatchPMapI;
typedef PatchPMap::const_iterator       PatchPMapCI;

typedef std::vector<PatchP>             PatchPVec;
typedef PatchPVec::iterator             PatchPVecI;
typedef PatchPVec::const_iterator       PatchPVecCI;

////////////////////////////////////////////////////////////////////////////////

/// Base class for patch objects.
///
///    A patch is a piece of 2D surface surrounding (part of) a 3D compartment.
///    This base class provides basic functionality and descriptive data that
///    is shared by all types of patches ('type' meaning different types of
///    geometric descriptions):
///
///        - Getting and setting a valid patch ID string, and handling
///          the interaction with the container object.
///
///        - Getting (and at least in this base class also setting) the total
///          area of the patch.
///
///        - The surface systems associated with the patches.
///
///        - References to inside/outside compartments.
///
///    This base class can be used directly with well-mixed solvers.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Patch
{
public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the patch.
    /// \param container Pointer to the parent geometry container.
    /// \param icomp Pointer to the inner compartment.
    /// \param ocomp Pointer to the outer compartment.
    /// \param surfsys Pointer to the associated surface system.
    /// \param area Area of the patch.
	Patch(std::string const & id, steps::wm::Geom * container,
		  steps::wm::Comp* icomp, steps::wm::Comp* ocomp = 0, double area = 0.0);

    /// Destructor
	virtual ~Patch(void);

	////////////////////////////////////////////////////////////////////////
	// PATCH PROPERTIES (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

	/// Return the patch id.
    ///
    /// \return ID of the patch.
	std::string getID(void) const
	{ return pID; }

	/// Set or change the patch id.
    ///
    /// \param id ID of the patch.
	void setID(std::string const & id);

	/// Return a pointer to the geometry container object.
    ///
    /// \return Pointer to the parent geometry container.
	steps::wm::Geom * getContainer(void) const
	{ return pContainer; }

    /// Return the area of the patch.
    ///
    /// \return Area of the patch.
	double getArea(void) const
	{ return pArea; }

    /// Set the area of the patch.
    ///
    /// \param area Area of the patch.
	virtual void setArea(double vol);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): VOLUME SYSTEM
	////////////////////////////////////////////////////////////////////////

    /// Add a surface system with name id.
    ///
    /// \param id ID of the surface system.
	void addSurfsys(std::string const & id);

    /// Get a surface system.
    ///
    /// \return List of the surface systems associated to the patch.
	std::set<std::string> getSurfsys(void) const
	{return pSurfsys; }

    /// Delete a surface system with name id.
    ///
    /// \param id ID of the surface system.
	void delSurfsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): COMPARTMENTS
	////////////////////////////////////////////////////////////////////////

	/// Return the inner compartment.
    ///
    /// \return Pointer to the inner compartment.
	steps::wm::Comp * getIComp(void) const
	{ return pIComp; }

    ///Return the outer compartment.
    ///
    /// \return Pointer to the outer compartment.
	steps::wm::Comp * getOComp(void) const
	{ return pOComp; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: PATCHES
	////////////////////////////////////////////////////////////////////////

	/// Set the inner compartment.
    ///
    /// \param icomp Pointer to the inner compartment.
	void _setIComp(steps::wm::Comp* icomp);

	/// Set the outer compartment.
    ///
    /// \param ocomp Pointer to the outer compartment.
	void _setOComp(steps::wm::Comp* ocomp);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
	/// Called if Python object deleted, or from del method in parent object.
	/// Will only be called once
	void _handleSelfDelete(void);

	////////////////////////////////////////////////////////////////////////

protected:

	////////////////////////////////////////////////////////////////////////

	double                              pArea;

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;
	steps::wm::Geom                   * pContainer;
	steps::wm::Comp                   * pIComp;
	steps::wm::Comp                   * pOComp;
	std::set<std::string>               pSurfsys;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wm)
END_NAMESPACE(steps)

#endif
// STEPS_WM_PATCH_HPP
