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

#ifndef STEPS_SOLVER_REACDEF_HPP
#define STEPS_SOLVER_REACDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/solver/api.hpp"
#include "steps/model/reac.hpp"
#include "steps/model/spec.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forwards declarations
class Statedef;

////////////////////////////////////////////////////////////////////////////////

/// Defined Reaction.
class Reacdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the reaction.
    /// \param r Pointer to the associated Reac object.
    Reacdef(Statedef * sd, uint idx, steps::model::Reac * r);

    /// Destructor
    ~Reacdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this reaction rule.
    inline uint gidx() const
    { return pIdx; }

    /// Return the name of the reaction.
    std::string const name() const;

    /// Return the order of this reaction.
    uint order() const;

    /// Return the MACROscopic reaction constant.
    double kcst() const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////
    /// \todo imcompleted.
    uint lhs(uint gidx) const;
    int dep(uint gidx) const;
    uint rhs(uint gidx) const;
    int upd(uint gidx) const;
    bool reqspec(uint gidx) const;

    inline steps::solver::gidxTVecCI bgnUpdColl() const
    { return pSpec_UPD_Coll.begin(); }
    inline steps::solver::gidxTVecCI endUpdColl() const
    { return pSpec_UPD_Coll.end(); }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;
    uint                                pIdx;
    std::string                         pName;
    uint                                 pOrder;
    double                                 pKcst;

    // The stoichiometry stored as model level Spec objects.
    // To be used during setup ONLY
    steps::model::SpecPVec                 pLhs;
    steps::model::SpecPVec                 pRhs;

    bool                                pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    int                               * pSpec_DEP;
    uint                              * pSpec_LHS;
    uint                              * pSpec_RHS;
    int                               * pSpec_UPD;
    steps::solver::gidxTVec             pSpec_UPD_Coll;

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_REACDEF_HPP

// END
