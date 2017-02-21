/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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

#ifndef STEPS_MODEL_CHAN_HPP
#define STEPS_MODEL_CHAN_HPP 1

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
class Chan;
class ChanState;

// Auxiliary declarations.
typedef Chan *                          ChanP;
typedef std::map<std::string, ChanP>    ChanPMap;
typedef ChanPMap::iterator              ChanPMapI;
typedef ChanPMap::const_iterator        ChanPMapCI;

typedef std::vector<ChanP>              ChanPVec;
typedef ChanPVec::iterator              ChanPVecI;
typedef ChanPVec::const_iterator        ChanPVecCI;

////////////////////////////////////////////////////////////////////////////////
/// Channel grouping a number of states with voltage-dependent
/// transitions permitted between conducting states.

///
/// \warning Methods start with an underscore are not exposed to Python.

class Chan
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the species.
    /// \param model Pointer to the parent model.
    Chan(std::string const & id, Model * model);

    /// Destructor
    ~Chan(void);

    ////////////////////////////////////////////////////////////////////////
    // CHANNEL PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the channel ID.
    ///
    /// \return ID of the channel.
    std::string getID(void) const
    { return pID; }

    /// Set or change the species ID.
    ///
    /// \param id ID of the species.
    void setID(std::string const & id);

    /// Return a pointer to the parent model.
    ///
    /// \return Pointer to the parent model.
    Model * getModel(void) const
    { return pModel; }

    ////////////////////////////////////////////////////////////////////////
    // OPERATIONS (EXPOSED TO PYTHON): CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    /// Get a channel state by its ID.
    ///
    /// \param id ID of the required channel state.
    /// \return Pointer to the channel state object.
    ChanState * getChanState(std::string const & id) const;

    /*
    /// Delete a channel state by its ID.
    ///
    /// \param id ID of the channel state to be deleted.
    void delChanState(std::string const & id);
    */

    /// Get all channel states stored in this channel.
    ///
    /// \return A vector of pointers to the channel state objects
    ///         stored in the system.
    std::vector<ChanState *> getAllChanStates(void) const;

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED) OPERATIONS: DELETION
    ////////////////////////////////////////////////////////////////////////

    /// Self delete.
    ///
    /// Called if Python object deleted, or from del method in parent object.
    /// Will only be called once
    void _handleSelfDelete(void);

    ////////////////////////////////////////////////////////////////////////
    // INTERNAL (NON-EXPOSED): STEPS::MODEL OPERATIONS
    ////////////////////////////////////////////////////////////////////////

    /// Add a state to the channel.
    ///
    /// \param cstate Pointer to the ChanState.
    void _handleChanStateAdd(ChanState * cstate);


    /// Delete a state from the channel.
    ///
    /// \param cstate Pointer to the ChanState.
    void _handleChanStateDel(ChanState * cstate);


    /// Check if a channel state id is occupied.
    ///
    /// \param id ID of the channel state.
    void _checkChanStateID(std::string const & id) const;

    /// Change the id of a reaction from o to n.
    ///
    /// \param o Old id of the channel state.
    /// \param n New id of the channel state.
    void _handleChanStateIDChange(std::string const & o, std::string const & n);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    std::string                              pID;
    Model                                  * pModel;
    std::map<std::string, ChanState *>       pChanStates;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MODEL_CHAN_HPP

// END
