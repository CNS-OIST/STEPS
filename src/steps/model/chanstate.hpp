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

#include "spec.hpp"

namespace steps::model {

/// Channel state.
/// Component that represents a channel state that can be referred to from
/// voltage-dependent transitions.
///
/// \warning Methods start with an underscore are not exposed to Python.

class ChanState: public Spec {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    /// Constructor
    ///
    /// \param id ID of the channel state.
    /// \param chan Reference to the parent channel.
    ChanState(std::string const& id, Model& model, Chan& chan);

    ChanState(const ChanState&) = delete;
    ChanState& operator=(const ChanState&) = delete;

    /// Destructor
    ~ChanState();

    ////////////////////////////////////////////////////////////////////////
    // CHANNEL STATE PROPERTIES
    ////////////////////////////////////////////////////////////////////////

    /// Return a reference to the associated channel.
    ///
    /// \return Reference to the associated channel.
    inline Chan& getChan() const noexcept {
        return pChan;
    }

    void setID(std::string const& id);

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

    Chan& pChan;
};

}  // namespace steps::model
