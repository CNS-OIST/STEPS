/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
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

#pragma once

#include <string>
#include <vector>

#include "fwd.hpp"
#include "model/fwd.hpp"

namespace steps::solver {

/// This class provides functionality for grouping a set of channel states
/// (which look just like Spec objects at this level). Channel states
/// behave like Spec objects that may be involved in voltage-dependent
/// transitions, and may diffuse in the membrane (in
/// which case the channel states must diffuse with it) dettach from the
/// membrane and diffuse in a volume (in which case channel state transitions
/// need to be turned off) for example. No need to involve Channels in this
/// description at all, it's all about the Channelstates that describe the
/// channel.

/// Defined Channel
class Chandef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the channel.
    /// \param c Reference to the associated Chan object.
    Chandef(Statedef& sd, chan_global_id idx, model::Chan& c);

    Chandef(const Chandef&) = delete;
    Chandef& operator=(const Chandef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this species.
    inline chan_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the species.
    inline const std::string& name() const noexcept {
        return pName;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    /// The global species indices of the associated channel states.
    inline const auto& chanstates() const noexcept {
        return pChanStates;
    }

    /// The number of channel states describing this channel.
    inline uint nchanstates() const noexcept {
        return pChanStatesVec.size();
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    void setup(const Statedef& sd);

  private:
    const chan_global_id pIdx;
    const std::string pName;
    bool pSetupdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    // The global indices of the channel states. Storing in arbitrary order
    // and storing only these indices rather than the usual big table to
    // see if I can get away with that.
    std::vector<spec_global_id> pChanStates;
    // Vector of the channel state objects
    // To be used during setup() ONLY
    const std::vector<model::ChanState*> pChanStatesVec;
};

}  // namespace steps::solver
