////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_MODEL_SURFSYS_HPP
#define STEPS_MODEL_SURFSYS_HPP 1


// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "../common.h"

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
/// Surface system.
/// Container that collects reactions involving a reactant
/// embedded in a membrane.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Surfsys
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the surface system.
    /// \param model Pointer to the parent model.
	Surfsys(std::string const & id, Model * model);

    /// Destructor
	~Surfsys(void);

	////////////////////////////////////////////////////////////////////////
	// SURFACE SYSTEM PROPERTIES
	////////////////////////////////////////////////////////////////////////

	/// Return the surface system ID.
    ///
    /// \return ID of the surface system.
	std::string getID(void) const
	{ return pID; }
	/// Set or change the surface system ID.
    ///
    /// \param id ID of the surface system.
	void setID(std::string const & id);

	/// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
	Model * getModel(void) const
	{ return pModel; }

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SURFACE REACTIONS
	////////////////////////////////////////////////////////////////////////

    /// Return a surface reaction with name id.
    ///
    /// \param id ID of the surface reaction.
    /// \return Pointer to the surface reaction.
	SReac * getSReac(std::string const & id) const;

    /// Delete a surace reaction with name id.
    ///
    /// \param id ID of the surface reaction.
	void delSReac(std::string const & id);

    /// Return a list of all surface reactions.
    ///
    /// \return List of pointers to surface reactions.
	std::vector<SReac *> getAllSReacs(void) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED) OPERATIONS: SURFACE REACTIONS
	////////////////////////////////////////////////////////////////////////

    /// Check if a surface reaction id is occupied.
    ///
    /// \param id ID of the surface reaction.
	void _checkSReacID(std::string const & id) const;

    /// Change a surface reaction id from o to n.
    ///
    /// \param o Old id of the surface reaction.
    /// \param n New id of the surface reaction.
	void _handleSReacIDChange(std::string const & o, std::string const & n);

    /// Add a surface reaction to the surface system.
    ///
    /// \param Pointer to the surface reaction.
	void _handleSReacAdd(SReac * sreac);

    /// Delete a surface reaction in the surface system.
    ///
    /// \param Pointer to the surface reaction.
	void _handleSReacDel(SReac * sreac);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS (EXPOSED TO PYTHON): SPECIES
	////////////////////////////////////////////////////////////////////////

    /// Return all species in the surface system.
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
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

    /// Count the surface reactions in the surface system.
    ///
    /// \return Number of surface reactions.
	inline uint _countSReacs(void) const
	{ return pSReacs.size(); }

    /// Get a surface reaction with index lidx.
    ///
    /// \param lidx Index of the surface reaction.
    /// \return Pointer to the surface reaction.
	SReac * _getSReac(uint lidx) const;

    /// Get all surface reactions in the surface system.
    ///
    /// \return Map of surface reactions.
	const std::map<std::string, SReac *> & _getAllSReacs(void) const
	{ return pSReacs; }

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
	////////////////////////////////////////////////////////////////////////
    /// Delete a species in the surface system.
    ///
    /// \param spec Pointer to the species.
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
