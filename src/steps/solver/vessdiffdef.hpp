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
#include <string>

// STEPS headers.
#include "model/vessdiff.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

/// Defined diffusion on vesicle surface object.

class VesSDiffdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the object.
    /// \param vsd Pointer to the associated VesSDiff object.
    VesSDiffdef(Statedef* sd, vessdiff_global_id idx, model::VesSDiff* vsd);

    /// Destructor
    ~VesSDiffdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this vesicle surface diffusion rule.
    inline vessdiff_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of this vesicle surface diffusion rule.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the diffusion constant.
    inline double dcst() const noexcept {
        return pDcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LIGAND
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of the ligand species.
    spec_global_id lig() const;

    int dep(spec_global_id gidx) const;
    bool reqspec(spec_global_id gidx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: DIFFUSION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Set the diffusion constant for this vesicle surface diffusion rule.
    ///
    /// \param d Rate constant of the vesicle surface diffusion rule.
    void setDcst(double d);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;

    // The global index of this diffusion rule
    vessdiff_global_id pIdx;

    // The string identifier of this diffusion rule
    std::string pName;

    // The diffusion constant
    double pDcst;

    // The chemical species to which this diffusion rule applies,
    // safer to store as a sting, rather than model level spec object pointer
    std::string pLig;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
