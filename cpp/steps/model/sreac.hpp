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

#ifndef STEPS_MODEL_SREAC_HPP
#define STEPS_MODEL_SREAC_HPP 1

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
class SReac;
class Surfsys;
class Model;
class Spec;

// Auxiliary declarations.
typedef SReac *						    SReacP;
typedef std::map<std::string, SReacP>   SReacPMap;
typedef SReacPMap::iterator             SReacPMapI;
typedef SReacPMap::const_iterator       SReacPMapCI;
typedef std::vector<SReacP>             SReacPVec;
typedef SReacPVec::iterator             SReacPVecI;
typedef SReacPVec::const_iterator       SReacPVecCI;

////////////////////////////////////////////////////////////////////////////////

class SReac
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

	SReac(std::string const & id, Surfsys * surfsys,
		  std::vector<Spec *> const & vlhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & slhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & irhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & srhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & orhs = std::vector<Spec *>(),
		  double kcst = 0.0);
	~SReac(void);

	////////////////////////////////////////////////////////////////////////
	// SURFACE REACTION RULE PROPERTIES
	////////////////////////////////////////////////////////////////////////

	// Return the surface reaction rule ID
	std::string getID(void) const
	{ return pID; }
	// Set or change the surface reaction rule ID
	void setID(std::string const & id);

	// Return a pointer to the parent surface system
	Surfsys * getSurfsys(void) const
	{ return pSurfsys; }

	// Return a pointer to the parent model
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON):
	////////////////////////////////////////////////////////////////////////

	bool getInner(void) const
	{ return (! pOuter); }
	void setInner(bool inner);

	bool getOuter(void) const
	{ return pOuter; }
	void setOuter(bool outer);

	const std::vector<Spec *> & getVLHS(void) const
	{ return pVLHS; }
	void setVLHS(std::vector<Spec *> const & vlhs);

	const std::vector<Spec *> & getSLHS(void) const
	{ return pSLHS; }
	void setSLHS(std::vector<Spec *> const & slhs);

	const std::vector<Spec *> & getIRHS(void) const
	{ return pIRHS; }
	void setIRHS(std::vector<Spec *> const & irhs);

	const std::vector<Spec *> & getSRHS(void) const
	{ return pSRHS; }
	void setSRHS(std::vector<Spec *> const & srhs);

	const std::vector<Spec *> & getORHS(void) const
	{ return pORHS; }
	void setORHS(std::vector<Spec *> const & orhs);

	uint getOrder(void) const
	{ return pOrder; }

	double getKcst(void) const
	{ return pKcst; }
	void setKcst(double kcst);

	// This method returns a list of all species involved in this
	// surface reaction, on both the left and righthand side and does not
	//contain any duplicate members
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
	Surfsys                           * pSurfsys;

	bool                                pOuter;
	std::vector<Spec *>                 pVLHS;
	std::vector<Spec *>                 pSLHS;
	std::vector<Spec *>                 pIRHS;
	std::vector<Spec *>                 pSRHS;
	std::vector<Spec *>                 pORHS;
	uint                                pOrder;
	double                              pKcst;

	////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(model)
END_NAMESPACE(steps)

#endif
// STEPS_MODEL_SREAC_HPP

// END
