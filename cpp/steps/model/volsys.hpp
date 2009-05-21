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

#ifndef STEPS_MODEL_VOLSYS_HPP
#define STEPS_MODEL_VOLSYS_HPP 1

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
class Volsys;
class Reac;
class Diff;

// Auxiliary declarations.
typedef Volsys *                        VolsysP;
typedef std::map<std::string, VolsysP>  VolsysPMap;
typedef VolsysPMap::iterator            VolsysPMapI;
typedef VolsysPMap::const_iterator      VolsysPMapCI;

////////////////////////////////////////////////////////////////////////////////

/// Container that collects reactions involving a reactant
/// enclosed in a compartment.
///
class Volsys
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Volsys(std::string const & id, Model * model);
	~Volsys(void);

	////////////////////////////////////////////////////////////////////////
	// VOLUME SYSTEM PROPERTIES
	////////////////////////////////////////////////////////////////////////

	/// Return the volume system ID.
	std::string getID(void) const
	{ return pID; }
	/// Set or change the volume system ID.
	void setID(std::string const & id);

	/// Return a pointer to the parent model.
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): REACTIONS
	////////////////////////////////////////////////////////////////////////

	Reac * getReac(std::string const & id) const;
	void delReac(std::string const & id);

	std::vector<Reac *> getAllReacs(void) const;

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): DIFFUSION
	////////////////////////////////////////////////////////////////////////

	Diff * getDiff(std::string const & id) const;
	void delDiff(std::string const & id);

	std::vector<Diff *> getAllDiffs(void) const;

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SPECIES
	////////////////////////////////////////////////////////////////////////

	// This method returns a list of all species involved in this
	// volume system, and does not contain any duplicate members
	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

	// Called if Python object deleted, or from del method in parent object.
	// Will only be called once
	void _handleSelfDelete(void);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: REACTIONS
	////////////////////////////////////////////////////////////////////////

	void _checkReacID(std::string const & id) const;
	void _handleReacIDChange(std::string const & o, std::string const & n);
	void _handleReacAdd(Reac * reac);
	void _handleReacDel(Reac * reac);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DIFFUSION
	////////////////////////////////////////////////////////////////////////

	void _checkDiffID(std::string const & id) const;
	void _handleDiffIDChange(std::string const & o, std::string const & n);
	void _handleDiffAdd(Diff * diff);
	void _handleDiffDel(Diff * diff);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

	inline uint _countReacs(void) const
	{ return pReacs.size(); }
	Reac * _getReac(uint lidx) const;

	const std::map<std::string, Reac *> & _getAllReacs(void) const
	{ return pReacs; }

	inline uint _countDiffs(void) const
	{ return pDiffs.size(); }
	Diff * _getDiff(uint lidx) const;

	const std::map<std::string, Diff *> & _getAllDiffs(void) const
	{ return pDiffs; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
	////////////////////////////////////////////////////////////////////////

	void _handleSpecDelete(Spec * spec);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;
	Model                             * pModel;
	std::map<std::string, Reac *>       pReacs;
	std::map<std::string, Diff *>       pDiffs;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_VOLSYS_HPP

// END
