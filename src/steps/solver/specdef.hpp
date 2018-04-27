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

#ifndef STEPS_SOLVER_SPECDEF_HPP
#define STEPS_SOLVER_SPECDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/model/spec.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forwards declarations
//class Statedef;

////////////////////////////////////////////////////////////////////////////////

/// Defined Species
class Specdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the species.
    /// \param d Pointer to the assocaited Spec object.
    Specdef(Statedef * sd, uint idx, steps::model::Spec * d);

    /// Destructor
    ~Specdef(void);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this species.
    inline uint gidx(void) const
    { return pIdx; }

    /// Return the name of the species.
    std::string const name(void) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    /// This method is included for consistency with other def objects,
    /// but currently does nothing.
    void setup(void);

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;
    uint                                pIdx;
    std::string                         pName;
    bool                                pSetupdone;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_SPECDEF_HPP

// END
