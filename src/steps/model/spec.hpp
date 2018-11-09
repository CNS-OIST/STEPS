/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
#    as published by the Free Software Foundation.
#    
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#    
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_MODEL_SPEC_HPP
#define STEPS_MODEL_SPEC_HPP 1

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace model {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Model;
class Spec;

// Auxiliary declarations.
typedef Spec *                          SpecP;
typedef std::map<std::string, SpecP>    SpecPMap;
typedef SpecPMap::iterator              SpecPMapI;
typedef SpecPMap::const_iterator        SpecPMapCI;

typedef std::vector<SpecP>              SpecPVec;
typedef SpecPVec::iterator              SpecPVecI;
typedef SpecPVec::const_iterator        SpecPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Species reactant.
/// Component that represents a reactant that can be referred to from
/// volume and surface systems.
///
/// \warning Methods start with an underscore are not exposed to Python.

class Spec
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the species.
    /// \param model Pointer to the parent model.
    Spec(std::string const & id, Model * model, int valence = 0);

    /// Destructor
    virtual ~Spec();

    ////////////////////////////////////////////////////////////////////////
    // SPECIES PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the species ID.
    ///
    /// \return ID of the species.
    std::string getID() const
    { return pID; }

    /// Set or change the species ID.
    ///
    /// \param id ID of the species.
    virtual void setID(std::string const & id);

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    Model * getModel() const
    { return pModel; }

    /// Set the valence of the species.
    ///
    /// \param valence Valence of the species.
    void setValence(int valence);

    /// Return the valence of the species.
    ///
    /// \return Valence of the species.
    int getValence() const
    { return pValence; }

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    virtual void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                         pID;
    Model                             * pModel;
    int                                 pValence;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_SPEC_HPP

// END
