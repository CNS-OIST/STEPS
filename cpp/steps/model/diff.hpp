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

#ifndef STEPS_MODEL_DIFF_HPP
#define STEPS_MODEL_DIFF_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <string>
#include <map>
#include <vector>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(model)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Diff;
class Volsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef Diff *						    DiffP;
typedef std::map<std::string, DiffP>    DiffPMap;
typedef DiffPMap::iterator              DiffPMapI;
typedef DiffPMap::const_iterator        DiffPMapCI;
typedef std::vector<DiffP>              DiffPVec;
typedef DiffPVec::iterator              DiffPVecI;
typedef DiffPVec::const_iterator        DiffPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Diff
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Diff(std::string const & id, Volsys * volsys, Spec * lig, double dcst=0.0);
	~Diff(void);

	////////////////////////////////////////////////////////////////////////
	// DIFFUSION RULE PROPERTIES
	////////////////////////////////////////////////////////////////////////

	// Return the diffusion rule ID.
	std::string getID(void) const
	{ return pID; }
	void setID(std::string const & id);

	// Return a pointer to the parent volume system.
	Volsys * getVolsys(void) const
	{ return pVolsys; }

	// Return a pointer to the parent model.
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON):
	////////////////////////////////////////////////////////////////////////

	// Return a pointer to the species to which this diffusion rule applies
	Spec * getLig(void) const
	{ return pLig; }
	void setLig(Spec * lig);

	double getDcst(void) const
	{ return pDcst; }
	void setDcst(double dcst);

	// Create a vector of all species in this diffusion rule
	// Currently will return only one species
	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////

	// Called if Python object deleted, or from del method in parent object.
	// Will only be called once
	void _handleSelfDelete(void);

	////////////////////////////////////////////////////////////////////////

private:

	////////////////////////////////////////////////////////////////////////

	std::string                         pID;
	Model                             * pModel;
	Volsys                            * pVolsys;
	Spec                              * pLig;
	double                              pDcst;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_DIFF_HPP

// END
