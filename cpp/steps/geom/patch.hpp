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
//
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_WM_PATCH_HPP
#define STEPS_WM_PATCH_HPP 1

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
#include <steps/model/surfsys.hpp>
#include <steps/geom/geom.hpp>
#include <steps/geom/comp.hpp>

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
///        * Getting and setting a valid patch ID string, and handling
///          the interaction with the container object.
///
///        * Getting (and at least in this base class also setting) the total
///          area of the patch.
///
///        * The surface systems associated with the patches.
///
///        * References to inside/outside compartments.
///
///    This base class can be used directly with well-mixed solvers.
///

class Patch
{
public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Patch(std::string const & id, steps::wm::Geom * container,
		  steps::wm::Comp* icomp, steps::wm::Comp* ocomp = 0,
		  steps::model::Surfsys* surfsys = 0, double area = 0.0);
	virtual ~Patch(void);

	////////////////////////////////////////////////////////////////////////
	// PATCH PROPERTIES (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

	// return the patch id
	std::string getID(void) const
	{ return pID; }
	// set or change the patch id
	void setID(std::string const & id);

	// return a pointer to the geometry container object
	steps::wm::Geom * getContainer(void) const
	{ return pContainer; }

	double getArea(void) const
	{ return pArea; }
	virtual void setArea(double vol);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): VOLUME SYSTEM
	////////////////////////////////////////////////////////////////////////

	void addSurfsys(std::string const & id);
	std::set<std::string> getSurfsys(void) const
	{return pSurfsys; }
	void delSurfsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// DATA ACCESS (EXPOSED TO PYTHON): COMPARTMENTS
	////////////////////////////////////////////////////////////////////////

	// return a reference to the inside compartment
	steps::wm::Comp * getIComp(void) const
	{ return pIComp; }
	steps::wm::Comp * getOComp(void) const
	{ return pOComp; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: PATCHES
	////////////////////////////////////////////////////////////////////////

	void _setIComp(steps::wm::Comp* icomp);
	void _setOComp(steps::wm::Comp* ocomp);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

	// Called if Python object deleted, or from del method in parent object.
	// Will only be called once
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
