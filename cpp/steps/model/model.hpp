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

#ifndef STEPS_MODEL_MODEL_HPP
#define STEPS_MODEL_MODEL_HPP 1

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
START_NAMESPACE(model)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Model;
class Spec;
class Surfsys;
class Volsys;
class Reac;
class SReac;
class Diff;

// Auxiliary declarations.

typedef Spec *                           SpecP;
typedef std::vector<SpecP>               SpecPVec;
typedef SpecPVec::iterator               SpecPVecI;
typedef SpecPVec::const_iterator         SpecPVecCI;

typedef Volsys *                        VolsysP;
typedef std::map<std::string, VolsysP>  VolsysPMap;
typedef VolsysPMap::iterator            VolsysPMapI;
typedef VolsysPMap::const_iterator      VolsysPMapCI;

typedef Surfsys *                       SurfsysP;
typedef std::map<std::string, SurfsysP> SurfsysPMap;
typedef SurfsysPMap::iterator           SurfsysPMapI;
typedef SurfsysPMap::const_iterator     SurfsysPMapCI;

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

/// Top-level container for the objects in a kinetic model.
///
/// A steps::model::Model object is parent to the following objects:
/// <UL>
/// <LI>steps::model::Spec
/// <LI>steps::model::Volsys
/// <LI>steps::model::Surfsys
/// </UL>
///

class Model
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Model(void);
	~Model(void);

	// Model * deepcopy(void);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: SPECIES (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

	Spec * getSpec(std::string const & id) const;
	void delSpec(std::string const & id);
	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: VOLSYS (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

	Volsys * getVolsys(std::string const & id) const;
	void delVolsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: SURFSYS (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

	Surfsys * getSurfsys(std::string const & id) const;
	void delSurfsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

	inline uint _countSpecs(void) const
	{ return pSpecs.size(); }
	Spec * _getSpec(uint gidx) const;

	uint _countReacs(void) const;
	Reac * _getReac(uint gidx) const;

	uint _countSReacs(void) const;
	SReac * _getSReac(uint gidx) const;

	uint _countDiffs(void) const;
	Diff * _getDiff(uint gidx) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
	////////////////////////////////////////////////////////////////////////

	void _checkSpecID(std::string const & id) const;
	void _handleSpecIDChange(std::string const & o, std::string const & n);
	void _handleSpecAdd(Spec * spec);
	void _handleSpecDel(Spec * spec);

	void _checkVolsysID(std::string const & id) const;
	void _handleVolsysIDChange(std::string const & o, std::string const & n);
	void _handleVolsysAdd(Volsys * volsys);
	void _handleVolsysDel(Volsys * volsys);

	void _checkSurfsysID(std::string const & id) const;
	void _handleSurfsysIDChange(std::string const & o, std::string const & n);
	void _handleSurfsysAdd(Surfsys * surfsys);
	void _handleSurfsysDel(Surfsys * surfsys);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::map<std::string, Spec *>       pSpecs;
	std::map<std::string, Volsys *>     pVolsys;
	std::map<std::string, Surfsys *>    pSurfsys;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_MODEL_HPP

// END
