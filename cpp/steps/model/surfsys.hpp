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

#ifndef STEPS_MODEL_SURFSYS_HPP
#define STEPS_MODEL_SURFSYS_HPP 1

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
class SReac;

// Auxiliary declarations.
typedef Surfsys *                       SurfsysP;
typedef std::map<std::string, SurfsysP> SurfsysPMap;
typedef SurfsysPMap::iterator           SurfsysPMapI;
typedef SurfsysPMap::const_iterator     SurfsysPMapCI;

////////////////////////////////////////////////////////////////////////////////

/// Container that collects reactions involving a reactant
/// embedded in a membrane.
///
class Surfsys
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Surfsys(std::string const & id, Model * model);
	~Surfsys(void);

	////////////////////////////////////////////////////////////////////////
	// SURFACE SYSTEM PROPERTIES
	////////////////////////////////////////////////////////////////////////

	/// Return the surface system ID.
	std::string getID(void) const
	{ return pID; }
	/// Set or change the surface system ID.
	void setID(std::string const & id);

	/// Return a pointer to the parent model.
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SURFACE REACTIONS
	////////////////////////////////////////////////////////////////////////

	SReac * getSReac(std::string const & id) const;
	void delSReac(std::string const & id);

	std::vector<SReac *> getAllSReacs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: SURFACE REACTIONS
	////////////////////////////////////////////////////////////////////////

	void _checkSReacID(std::string const & id) const;
	void _handleSReacIDChange(std::string const & o, std::string const & n);
	void _handleSReacAdd(SReac * sreac);
	void _handleSReacDel(SReac * sreac);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SPECIES
	////////////////////////////////////////////////////////////////////////

	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

	// Called if Python object deleted, or from del method in parent object.
	// Will only be called once
	void _handleSelfDelete(void);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

	inline uint _countSReacs(void) const
	{ return pSReacs.size(); }
	SReac * _getSReac(uint lidx) const;

	const std::map<std::string, SReac *> & _getAllSReacs(void) const
	{ return pSReacs; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
	////////////////////////////////////////////////////////////////////////

	void _handleSpecDelete(Spec * spec);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;
	Model                             * pModel;
	std::map<std::string, SReac *>      pSReacs;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_SURFSYS_HPP

// END
