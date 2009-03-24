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

#ifndef STEPS_MODEL_REAC_HPP
#define STEPS_MODEL_REAC_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <cassert>
#include <string>
#include <vector>
#include <map>

// STEPS headers.
#include <steps/common.h>

////////////////////////////////////////////////////////////////////////////////

START_NAMESPACE(steps)
START_NAMESPACE(model)

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Reac;
class Volsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef Reac *						     ReacP;
typedef std::map<std::string, ReacP>     ReacPMap;
typedef ReacPMap::iterator               ReacPMapI;
typedef ReacPMap::const_iterator         ReacPMapCI;
typedef std::vector<ReacP>               ReacPVec;
typedef ReacPVec::iterator               ReacPVecI;
typedef ReacPVec::const_iterator         ReacPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Reac
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	Reac(std::string const & id, Volsys * volsys,
		 std::vector<Spec *> const & lhs = std::vector<Spec *>(),
		 std::vector<Spec *> const & rhs = std::vector<Spec *>(),
		 double kcst = 0.0);
	~Reac(void);

	////////////////////////////////////////////////////////////////////////
	// REACTION RULE PROPERTIES
	////////////////////////////////////////////////////////////////////////

	// Return the reaction rule ID.
	std::string getID(void) const
	{ return pID; }
	// Set or change the reaction rule ID.
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

	std::vector<Spec *> getLHS(void) const
	{ return pLHS; }
	void setLHS(std::vector<Spec *> const & lhs);

	std::vector<Spec *> getRHS(void) const
	{ return pRHS; }
	void setRHS(std::vector<Spec *> const & rhs);

	// This method returns a list of all species involved in this
	// reaction, on both the left and right-hand side and does not contain
	// any duplicate members
	std::vector<Spec *> getAllSpecs(void) const;

	uint getOrder(void) const
	{ return pOrder; }

	double getKcst(void) const
	{ return pKcst; }

	void setKcst(double kcst);

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

	std::vector<Spec *>                 pLHS;
	std::vector<Spec *>                 pRHS;
	uint                                pOrder;
	double                              pKcst;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_REAC_HPP

// END
