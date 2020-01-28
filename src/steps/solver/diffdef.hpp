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


#ifndef STEPS_SOLVER_DIFFDEF_HPP
#define STEPS_SOLVER_DIFFDEF_HPP 1


// STL headers.
#include <string>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/statedef.hpp"
#include "steps/model/diff.hpp"

////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace solver {

// Forwards declarations
//class Statedef;

////////////////////////////////////////////////////////////////////////////////
/// Defined diffusion object.
class Diffdef
{

public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param d Pointer to the associated Diff object.
    Diffdef(Statedef * sd, uint idx, steps::model::Diff * d);

    /// Destructor
    ~Diffdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this diffusion rule.
    inline uint gidx() const noexcept
    { return pIdx; }

    /// Return the name of this diffusion rule.
    inline const std::string& name() const noexcept
    { return pName; }

    /// Return the diffusion constant.
    inline double dcst() const noexcept
    { return pDcst; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LIGAND
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of the ligand species.
    uint lig() const;

    /// \todo check
    int dep(uint gidx) const;
    /// \todo check
    bool reqspec(uint gidx) const;


    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Set the diffusion constant for this diffusion rule.
    ///
    /// \param d Rate constant of the diffusion rule.
    void setDcst(double d);

    /*
    /// Set the ligand of the diffusion rule by its global index
    ///
    /// \param gidx Global index of the diffusion rule.
    void setLig(uint gidx);
    */

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    Statedef                          * pStatedef;

    // The global index of this diffusion rule
    uint                                pIdx;

    // The string identifier of this diffusion rule
    std::string                         pName;

    // The diffusion constant
    double                              pDcst;

    // The chemical species to which this diffusion rule applies,
    // safer to store as a sting, rather than model level spec object pointer
    std::string                         pLig;

    // The global index of the spec
    uint                                ligGIdx;

    bool                                pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: LIGAND
    ////////////////////////////////////////////////////////////////////////

    int                               * pSpec_DEP;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_SOLVER_DIFFDEF_HPP

// END

