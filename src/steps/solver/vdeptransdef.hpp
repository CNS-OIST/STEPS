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


#ifndef STEPS_SOLVER_VDEPTRANSDEF_HPP
#define STEPS_SOLVER_VDEPTRANSDEF_HPP 1

// STL headers.
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/model/vdeptrans.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forward declarations.
class VDepTransdef;

// Auxiliary declarations.
typedef VDepTransdef *                   VDepTransDefP;
typedef std::vector<VDepTransDefP>       VDepTransDefPVec;
typedef VDepTransDefPVec::iterator       VDepTransDefPVecI;
typedef VDepTransDefPVec::const_iterator VDepTransDefPVecCI;

////////////////////////////////////////////////////////////////////////////////

class VDepTransdef
{

public:

    /// Constructor
    ///
    /// \param sd Defined state of the solver.
    /// \param idx Global index of the voltage-dependent transition.
    /// \param vdt Pointer to the VDepTrans object.
    VDepTransdef(Statedef * sd, uint gidx, steps::model::VDepTrans * vdt);

    /// Destructor
    ~VDepTransdef();

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
    /*
    void setSrcChan(uint gidx);
    void setDstChan(uint gidx);
     */
    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT TRANSITION
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this voltage-dependent transition.
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of the voltage-dependent transition.
    inline std::string const name() const noexcept
    { return pName; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDANCE
    ////////////////////////////////////////////////////////////////////////

    /// Returns the transition rate for value of V in the range.
    ///
    double getVDepRate(double v) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of the source channel state (channel states
    /// stored as species objects)
    uint srcchanstate() const;

    /// Return the global index of the destination channel
    uint dstchanstate() const;

    int dep(uint gidx) const;
    bool req(uint gidx) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;

    // The global index of this voltage-dependent transition.
    uint                                pIdx;

    // The string identifier of this voltage-dependent transition.
    std::string                         pName;

    // True if setup() has been called.
    bool                                pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VOLTAGE DEPENDENCE
    ////////////////////////////////////////////////////////////////////////

    // The minimum voltage of stored voltage-dependent rates
    double                              pVMin;

    // The maximum voltage of stored voltage-dependent rates
    double                              pVMax;

    // The step between stored voltage-dependent rates
    double                              pDV;

    // Table of voltage-dependent rates, size (pVMax-pVMin)/pDV
    double                               * pVRateTab;

    ////////////////////////////////////////////////////////////////////////
    // DATA: CHANNEL STATES
    ////////////////////////////////////////////////////////////////////////

    // The channel state describing the 'source' channel
    // safer to store as a string, rather than model level spec object pointer
    std::string                         pSrc;

    // The channel state describing the 'destination' channel
    // safer to store as a string, rather than model level spec object pointer
    std::string                         pDst;

    int                               * pSpec_DEP;

    // Global index of the 'source' channel state
    uint                                pSpec_SRCCHAN;

    // Global index of the 'destination' channel state
    uint                                pSpec_DSTCHAN;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_VDEPTRANSDEF_HPP

// END
