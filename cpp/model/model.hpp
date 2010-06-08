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

#ifndef STEPS_MODEL_MODEL_HPP
#define STEPS_MODEL_MODEL_HPP 1


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
/// <LI>steps::model::Diff
/// </UL>
/// \sa Spec, Volsys, Surfsys, Diff.
/// \warning Methods start with an underscore are not exposed to Python.
///

class Model
{

public:

	////////////////////////////////////////////////////////////////////////
	// OBJECT CONSTRUCTION & DESTRUCTION
	////////////////////////////////////////////////////////////////////////

    /// Constructor
	Model(void);
    /// Destructor
	~Model(void);

	// Model * deepcopy(void);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: SPECIES (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////
    /// Return a species with name id.
    ///
    /// \param id ID of the species.
    /// \return Pointer to the species.
	Spec * getSpec(std::string const & id) const;

    /// Delete a species with name id.
    ///
    /// \param id ID of the species.
	void delSpec(std::string const & id);

    /// Return a list of all species in the Model object.
    ///
    /// \return List of pointers to the species in the Model object.
	std::vector<Spec *> getAllSpecs(void) const;

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: VOLSYS (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

    /// Return a volume system with name id.
    ///
    /// \param id ID of the volume system.
    /// \return Pointer to the volume system.
	Volsys * getVolsys(std::string const & id) const;

    /// Delete a volume system with name id.
    ///
    /// \param id ID of the volume system.
	void delVolsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// OPERATIONS: SURFSYS (EXPOSED TO PYTHON)
	////////////////////////////////////////////////////////////////////////

    /// Return a surface system with name id.
    ///
    /// \param id ID of the surface system.
    /// \return Pointer to the surface system.
	Surfsys * getSurfsys(std::string const & id) const;

    /// Delete a surface system with name id.
    ///
    /// \param id ID of the surface system.
	void delSurfsys(std::string const & id);

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
	////////////////////////////////////////////////////////////////////////

    /// Count the species in the Model object.
    ///
    /// \return Number of species.
	inline uint _countSpecs(void) const
	{ return pSpecs.size(); }

    /// Return a species with index gidx.
    ///
    /// \param gidx Index of the species.
    /// \return Pointer to the species.
	Spec * _getSpec(uint gidx) const;

    /// Count the reactions in the Model object.
    ///
    /// \return Number of reactions.
	uint _countReacs(void) const;

    /// Return a reaction with index gidx
    ///
    /// \param gidx Index of the reaction.
    /// \return Pointer to the reaction.
	Reac * _getReac(uint gidx) const;

    /// Count the surface reactions in the Model object.
    ///
    /// \return Number of surface reactions.
	uint _countSReacs(void) const;

    /// Return a surface with index gidx.
    ///
    /// \param gidx Index of the surface reaction.
    /// \return Pointer to the surface reaction.
	SReac * _getSReac(uint gidx) const;

    /// Count the diffusions in the Model object.
    ///
    /// \return Number of diffusions.
	uint _countDiffs(void) const;

    /// Return a diffusion with index gidx.
    ///
    /// \param gidx Index of the diffusion.
    /// \return Pointer to the diffusion.
	Diff * _getDiff(uint gidx) const;

	////////////////////////////////////////////////////////////////////////
	// INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
	////////////////////////////////////////////////////////////////////////

    /// Check if a species id is occupied.
    ///
    /// \param id ID of the species.
	void _checkSpecID(std::string const & id) const;

    /// Change the id of a species from o to n.
    ///
    /// \param o Old id of the species.
    /// \param n New id of the species.
	void _handleSpecIDChange(std::string const & o, std::string const & n);

    /// Add a species to the Model.
    ///
    /// \param spec Pointer to the species being added.
	void _handleSpecAdd(Spec * spec);

    /// Delete a species in the Model.
    ///
    /// \param spec Pointer to the species being deleted.
	void _handleSpecDel(Spec * spec);

    /// Check if a volume system id is occupied.
    ///
    /// \param id ID of the volume system.
	void _checkVolsysID(std::string const & id) const;

    /// Change the id of a volume system from o to n.
    ///
    /// \param o Old id of the volume system.
    /// \param n New id of the volume system.
	void _handleVolsysIDChange(std::string const & o, std::string const & n);

    /// Add a volume system to the Model.
    ///
    /// \param volsys Pointer to the volume system being added.
	void _handleVolsysAdd(Volsys * volsys);

    /// Delete a volume system in the Model.
    ///
    /// \param volsys Pointer to the volume system being deleted.
    void _handleVolsysDel(Volsys * volsys);

    /// Check if a surface system id is occupied.
    ///
    /// \param id ID of the surface system.
	void _checkSurfsysID(std::string const & id) const;

    /// Change the id of a surface system from o to n.
    ///
    /// \param o Old id of the surface system.
    /// \param n New id of the surface system.
	void _handleSurfsysIDChange(std::string const & o, std::string const & n);

    /// Add a surface system to the Model.
    ///
    /// \param surfsys Pointer to the surface system being added.
	void _handleSurfsysAdd(Surfsys * surfsys);

    /// Delete a surface system in the Model.
    ///
    /// \param surfsys Pointer to the surface system being deleted.
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
