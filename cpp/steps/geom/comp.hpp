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

#ifndef STEPS_WM_COMP_HPP
#define STEPS_WM_COMP_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>
#include <set>

// STEPS headers.
#include <steps/common.h>
#include <steps/model/volsys.hpp>

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
///        * Getting and setting a valid compartment ID string, and handling
///          the interaction with the container object.
///        * Getting (and at least in this base class also setting) the total
///          volume of the compartment.
///        * The volume systems implemented in the compartment.
///        * References to Patch objects.
///
///    This base class can be used directly with well-mixed solvers.
///

class Comp
{
public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Comp(std::string const & id, steps::wm::Geom * container,
		 steps::model::Volsys * volsys = 0, double vol = 0.0);
	virtual ~Comp(void);

	////////////////////////////////////////////////////////////////////////
	// COMPARTMENT PROPERTIES
	////////////////////////////////////////////////////////////////////////

	// return the compartment id
	std::string getID(void) const
	{ return pID; }
	// set or change the compartment id
	void setID(std::string const & id);

	// return a pointer to the geometry container object
	steps::wm::Geom * getContainer(void) const
	{ return pContainer; }

	double getVol(void) const
	{ return pVol; }
	virtual void setVol(double vol);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): VOLUME SYSTEM
	////////////////////////////////////////////////////////////////////////

	void addVolsys(std::string const & id);
	std::set<std::string> getVolsys(void) const
	{ return pVolsys; }
	void delVolsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): PATCHES
	////////////////////////////////////////////////////////////////////////

	// return a copy of the set of internal patches
	std::set<steps::wm::Patch *> getIPatches(void) const
	{ return pIPatches; }
	// return a copy of the set of external patches
	std::set<steps::wm::Patch *> getOPatches(void) const
	{ return pOPatches; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: PATCHES
	////////////////////////////////////////////////////////////////////////

	void _addIPatch(steps::wm::Patch * patch);
	void _delIPatch(steps::wm::Patch * patch);

	void _addOPatch(steps::wm::Patch * patch);
	void _delOPatch(steps::wm::Patch * patch);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

	// Called if Python object deleted, or from del method in parent object.
	// Will only be called once
	void _handleSelfDelete(void);


	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

	inline uint _countIPatches(void) const
	{ return pIPatches.size(); }

	steps::wm::Patch * _getIPatch(uint lidx) const;

	inline uint _countOPatches(void) const
	{ return pIPatches.size(); }

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
