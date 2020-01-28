/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_MODEL_CHANSTATE_HPP
#define STEPS_MODEL_CHANSTATE_HPP 1

// STL headers.
#include <cassert>
#include <map>
#include <string>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/model/spec.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace model {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations.
class Chan;
class ChanState;

// Auxiliary declarations.
typedef ChanState *                          ChanStateP;
typedef std::map<std::string, ChanStateP>    ChanStatePMap;
typedef ChanStatePMap::iterator              ChanStatePMapI;
typedef ChanStatePMap::const_iterator        ChanStatePMapCI;

typedef std::vector<ChanStateP>              ChanStatePVec;
typedef ChanStatePVec::iterator              ChanStatePVecI;
typedef ChanStatePVec::const_iterator        ChanStatePVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Channel state.
/// Component that represents a channel state that can be referred to from
/// voltage-dependent transitions.
///
/// \warning Methods start with an underscore are not exposed to Python.

class ChanState: public Spec
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the channel state.
    /// \param chan Pointer to the parent channel.
    ChanState(std::string const & id, Model* model, Chan * chan);

    /// Destructor
    ~ChanState();

    ////////////////////////////////////////////////////////////////////////
    // CHANNEL STATE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return a pointer to the associated channel.
    ///
    /// \return Pointer to the associated channel.
    inline Chan * getChan() const noexcept
    { return pChan; }

    void setID(std::string const & id);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete();

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): SOLVER HELPER METHODS
    ////////////////////////////////////////////////////////////////////////

    // ...

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Chan                              * pChan;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_CHANSTATE_HPP

// END

