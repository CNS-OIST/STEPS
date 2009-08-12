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
/// Surface reaction.
///
/// A SReac object describes a reaction which takes place on a surface system, 
/// i.e. a patch between two compartments.
///
/// \warning Methods start with an underscore are not exposed to Python.
class SReac
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////
    /// Constructor
    ///
    /// \param id ID of the surface reaction.
    /// \param surfsys Pointer to the parent surface system.
    /// \param vlhs Volume species on the left hand side of the reaction.
    /// \param slhs Surface species on the left hand side of the reaction.
    /// \param irhs Volume species in the inner compartment
    ///             and on the right hand side of the reaction.
    /// \param srhs Surface species on the right hand side of the reaction.
    /// \param orhs Volume species in the outer compartment
    ///             and on the right hand side of the reaction.
    /// \param kcst Rate constant of the reaction.
    ///
    /// \warning By default, the vlhs are defined in the outer compartment.
    ///          call setInner and SetOuter to change this default setting.
    /// \sa setInner, setOuter.
	SReac(std::string const & id, Surfsys * surfsys,
		  std::vector<Spec *> const & vlhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & slhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & irhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & srhs = std::vector<Spec *>(),
		  std::vector<Spec *> const & orhs = std::vector<Spec *>(),
		  double kcst = 0.0);

    /// Destructor
	~SReac(void);

	////////////////////////////////////////////////////////////////////////
	// SURFACE REACTION RULE PROPERTIES
	////////////////////////////////////////////////////////////////////////

	/// Return the surface reaction rule ID.
    ///
    /// \return ID of the surface reaction.
	std::string getID(void) const
	{ return pID; }

	/// Set or change the surface reaction rule ID.
    ///
    /// \param id ID of the surface reaction.
	void setID(std::string const & id);

	/// Return a pointer to the parent surface system.
    ///
    /// \return Pointer to the surface system.
	Surfsys * getSurfsys(void) const
	{ return pSurfsys; }

	/// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON):
	////////////////////////////////////////////////////////////////////////

    /// Check if the vlhs is assigned to the inner compartment.
    ///
    /// \return True if vlhs is assigned to the inner compartment;
    ///         False if else.   
	bool getInner(void) const
	{ return (! pOuter); }

    /// Assign the vlhs to the inner compartment.
    ///
    /// \param inner True to assign vlhs to the inner compartment;
    ///              False to assign vlhs to the outer compartment.
	void setInner(bool inner);

    /// Check if the vlhs is assigned to the outer compartment.
    ///
    /// \return True if vlhs is assigned to the outer compartment;
    ///         False if else.
	bool getOuter(void) const
	{ return pOuter; }

    /// Assign the vlhs to the outer compartment.
    ///
    /// \param outer True to assign vlhs to the outer compartment;
    ///              False to assign vlhs to the outer compartment.
	void setOuter(bool outer);

    /// Return a list of volume species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side volume species.
	const std::vector<Spec *> & getVLHS(void) const
	{ return pVLHS; }

    /// Set the volume species on the left hand side of reaction.
    ///
    /// \param vlhs Volume species on the left hand side of reaction.
	void setVLHS(std::vector<Spec *> const & vlhs);

    /// Return a list of surface species on the left hand side of reaction.
    ///
    /// \return List of pointers of left hand side surface species.
	const std::vector<Spec *> & getSLHS(void) const
	{ return pSLHS; }

    /// Set the surface species on the left hand side of reaction.
    ///
    /// \param slhs Surface species on the left hand side of reaction.
	void setSLHS(std::vector<Spec *> const & slhs);

    /// Return a list of inner volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side inner volume species.
	const std::vector<Spec *> & getIRHS(void) const
	{ return pIRHS; }

    /// Set the inner volume species on the right hand side of reaction.
    ///
    /// \param irhs Inner volume species on the right hand side of reaction.
	void setIRHS(std::vector<Spec *> const & irhs);

    /// Return a list of surface species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side surface species.
	const std::vector<Spec *> & getSRHS(void) const
	{ return pSRHS; }

    /// Set the surface species on the right hand side of reaction.
    ///
    /// \param srhs Surface species on the right hand side of reaction.
	void setSRHS(std::vector<Spec *> const & srhs);

    /// Return a list of outer volume species on the right hand side of reaction.
    ///
    /// \return List of pointers of right hand side outer volume species.
	const std::vector<Spec *> & getORHS(void) const
	{ return pORHS; }

    /// Set the outer volume species on the right hand side of reaction.
    ///
    /// \param orhs Outer volume species on the right hand side of reaction.
	void setORHS(std::vector<Spec *> const & orhs);

    /// Get the order of the surface reaction.
    ///
    /// \return Order of the reaction.
	uint getOrder(void) const
	{ return pOrder; }

    /// Get the rate constant of the surface reaction.
    ///
    /// \return Rate constant of the surface reaction.
	double getKcst(void) const
	{ return pKcst; }

    /// Set the rate constant of the surface reaction.
    ///
    /// \param Rate constant of the surface reaction.
	void setKcst(double kcst);

    /// Get a list of all species.
	///
    /// Returns a list of all species involved in this
	/// surface reaction, on both the left and righthand side
    /// and does not contain any duplicate members.
    ///
    /// \return List of pointers to the species.
	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
	////////////////////////////////////////////////////////////////////////
    /// Self delete.
    ///
	/// Called if Python object deleted, or from del method in parent object.
	/// Will only be called once
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
