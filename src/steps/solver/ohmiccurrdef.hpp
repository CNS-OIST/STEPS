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

// STL headers.
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// STEPS headers.
#include "model/ohmiccurr.hpp"
#include "statedef.hpp"
#include "types.hpp"
#include "util/common.hpp"

namespace steps::solver {

// Forward declarations.
class OhmicCurrdef;

// Auxiliary declarations.
typedef OhmicCurrdef* OhmicCurrDefP;
typedef std::vector<OhmicCurrDefP> OhmicCurrDefPVec;
typedef OhmicCurrDefPVec::iterator OhmicCurrDefPVecI;
typedef OhmicCurrDefPVec::const_iterator OhmicCurrDefPVecCI;

class OhmicCurrdef {
  public:
    /// Constructor
    ///
    /// \param sd Defined state of the solver.
    /// \param gidx Global index of the ohmic current.
    /// \param oc Pointer to the OhmicCurr object.
    OhmicCurrdef(Statedef* sd, ohmiccurr_global_id gidx, model::OhmicCurr* oc);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this ohmic current.
    inline ohmiccurr_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the ohmic current.
    inline const std::string& name() const noexcept {
        return pName;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PARAMETERS
    ////////////////////////////////////////////////////////////////////////

    inline double getG() const noexcept {
        return pG;
    }
    inline double getERev() const noexcept {
        return pERev;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL STATE
    ////////////////////////////////////////////////////////////////////////

    // Return the global index of the channel state
    spec_global_id chanstate() const;

    int dep(spec_global_id gidx) const;
    bool req(spec_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;

    // The global index of this ohmic current.
    ohmiccurr_global_id pIdx;

    // The string identifier of this ohmic current.
    std::string pName;

    // True if setup() has been called.
    bool pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: PARAMETERS
    ////////////////////////////////////////////////////////////////////////

    // The channel state store as a string, rather than ChanState pointer
    std::string pChanState;

    // Single-channel conductance
    double pG;

    // Reversal potential
    double pERev;

    util::strongid_vector<spec_global_id, int> pSpec_DEP;

    // Global index of the channel state
    spec_global_id pSpec_CHANSTATE;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
