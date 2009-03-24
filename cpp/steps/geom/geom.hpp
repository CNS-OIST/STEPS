////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WM_GEOM_HPP
#define STEPS_WM_GEOM_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(wm)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Comp;
class Patch;

// Auxiliary declarations.


////////////////////////////////////////////////////////////////////////////////

/// Test whether id is a valid identifier for some named STEPS component.
/// It returns a boolean value.
///
/// A valid id must be at least 1 character long and must start with an
/// underscore or an alphabetical character (a to z, or A to Z). Following
/// characters can be any alphanumerical character or again underscores.
/// The id cannot contain spaces.
///
/// Examples of valid id's:
///     a, _a, _a_, a000, adasf0, FSDaa9
///
STEPS_EXTERN bool isValidID(std::string const & id);

/// Test whether id is a valid identifier for some named STEPS component.
///
/// This function calls steps::model::isValidID() for the test. But whereas
/// isValidID() returns true or false, checkID() raises a steps::ArgErr
/// exception if the id is not valid. This makes it useful for checking
/// input arguments.
///
STEPS_EXTERN void checkID(std::string const & id);

////////////////////////////////////////////////////////////////////////////////

class Geom
{
public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Geom(void);
	virtual ~Geom(void);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON) : COMPARTMENTS
	////////////////////////////////////////////////////////////////////////

	steps::wm::Comp * getComp(std::string const & id) const;
	void delComp(std::string const & id);
	std::vector<steps::wm::Comp *> getAllComps(void) const;

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): PATCHES
	////////////////////////////////////////////////////////////////////////

	steps::wm::Patch * getPatch(std::string const & id) const;
	void delPatch(std::string const & id);
	std::vector<steps::wm::Patch *> getAllPatches(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

	inline uint _countComps(void) const
	{ return pComps.size(); }

	steps::wm::Comp * _getComp(uint gidx) const;

	inline uint _countPatches(void) const
	{ return pPatches.size(); }

	steps::wm::Patch * _getPatch(uint gidx) const;


	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::WM OPERATIONS
	////////////////////////////////////////////////////////////////////////

    void _checkCompID(std::string const & id) const;
    void _handleCompIDChange(std::string const & o, std::string const & n);
    void _handleCompAdd(steps::wm::Comp * comp);
    void _handleCompDel(steps::wm::Comp * comp);

    void _checkPatchID(std::string const & id) const;
    void _handlePatchIDChange(std::string const & o, std::string const & n);
    void _handlePatchAdd(steps::wm::Patch * patch);
    void _handlePatchDel(steps::wm::Patch *patch);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::map<std::string, steps::wm::Comp *>       pComps;
	std::map<std::string, steps::wm::Patch *>      pPatches;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(wm)
END_NAMESPACE(steps)

#endif
// STEPS_WM_GEOM_HPP

// END
