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


#ifndef STEPS_SOLVER_OHMICCURRDEF_HPP
#define STEPS_SOLVER_OHMICCURRDEF_HPP 1

// STL headers.
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/model/ohmiccurr.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forward declarations.
class OhmicCurrdef;

// Auxiliary declarations.
typedef OhmicCurrdef *                   OhmicCurrDefP;
typedef std::vector<OhmicCurrDefP>       OhmicCurrDefPVec;
typedef OhmicCurrDefPVec::iterator       OhmicCurrDefPVecI;
typedef OhmicCurrDefPVec::const_iterator OhmicCurrDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class OhmicCurrdef
{

public:
    /// Constructor
    ///
    /// \param sd Defined state of the solver.
    /// \param gidx Global index of the ohmic current.
    /// \param oc Pointer to the OhmicCurr object.
    OhmicCurrdef(Statedef * sd, uint gidx, steps::model::OhmicCurr * oc);

    /// Destructor
    ~OhmicCurrdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);


    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this ohmic current.
    inline uint gidx() const
    { return pIdx; }

    /// Return the name of the ohmic current.
    inline std::string const name() const
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PARAMETERS
    ////////////////////////////////////////////////////////////////////////

    inline double getG() const
    { return pG; }
    inline double getERev() const
    { return pERev; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL STATE
    ////////////////////////////////////////////////////////////////////////

    // Return the global index of the channel state
    uint chanstate() const;

    int dep(uint gidx) const;
    bool req(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;

    // The global index of this ohmic current.
    uint                                pIdx;

    // The string identifier of this ohmic current.
    std::string                         pName;

    // True if setup() has been called.
    bool                                pSetupdone;


    ////////////////////////////////////////////////////////////////////////
    // DATA: PARAMETERS
    ////////////////////////////////////////////////////////////////////////

    // The channel state store as a string, rather than ChanState pointer
    std::string                         pChanState;

    // Single-channel conductance
    double                              pG;

    // Reversal potential
    double                              pERev;


    int *                                  pSpec_DEP;

    // Global index of the channel state
    uint                                   pSpec_CHANSTATE;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_OHMICCURRDEF_HPP

// END
