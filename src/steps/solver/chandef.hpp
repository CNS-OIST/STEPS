/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_SOLVER_CHANDEF_HPP
#define STEPS_SOLVER_CHANDEF_HPP 1

// STL headers.
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/model/chan.hpp"
#include "steps/model/chanstate.hpp"


////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forwards declarations

////////////////////////////////////////////////////////////////////////////////

// TODO: Figure out what data this class will need to hold about child
// ChanStates. It is intended that the simulation solver will access this object
// meaning it will need to compute the populations of the ChannelStates depending
// on the interactions between them, but 1st version may not include this feature.
// In that case do I need this object?

/// This class provides functionality for grouping a set of channel states
/// (which look just like Spec objects at this level). Channel states
/// behave like Spec objects that may be involved in voltage-dependent
/// transitions, and may diffuse in the membrane (in
/// which case the channel states must diffuse with it) dettach from the membrane
/// and diffuse in a volume (in which case channel state transitions need to be
/// turned off) for example. No need to involve Channels in this description at all,
/// it's all about the Channelstates that describe the channel.

/// Defined Channel
class Chandef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the channel.
    /// \param c Pointer to the associated Chan object.
    Chandef(Statedef * sd, uint idx, steps::model::Chan * c);

    /// Destructor
    ~Chandef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this species.
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of the species.
    inline const std::string& name() const noexcept
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    // The global species indices of the associated channel states.
    inline uint * chanstates() const noexcept
    { return pChanStates; }

    // The number of channel states describing this channel.
    inline uint nchanstates() const noexcept
    { return pNChanStates; }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    void setup();

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;
    uint                                pIdx;
    std::string                         pName;
    bool                                pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    // The global indices of the channel states. Storing in arbitrary order
    // and storing only these indices rather than the usual big table to
    // see if I can get away with that.
    uint                              * pChanStates;
    uint                                 pNChanStates;
    // Vector of the channel state objects
    // To be used during setup() ONLY
    steps::model::ChanStatePVec            pChanStatesVec;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif

// STEPS_SOLVER_CHANDEF_HPP

// END


