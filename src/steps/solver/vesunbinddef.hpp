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
#include "model/spec.hpp"
#include "model/vesicle.hpp"
#include "model/vesunbind.hpp"
#include "solver/api.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

/// Defined vesicle-vesicle unbinding reaction.
class VesUnbinddef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the vesicle binding reaction.
    /// \param vb Pointer to the associated VesUnbind object.
    VesUnbinddef(Statedef* sd, vesunbind_global_id idx, model::VesUnbind* vb);

    /// Destructor
    ~VesUnbinddef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESUNBINDTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this vesicle binding reaction rule.
    inline vesunbind_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the vesicle binding reaction.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the MACROscopic vesicle binding reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    inline int immobility() const noexcept {
        return pImmobility;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    vesicle_global_id getVes1idx() const;
    vesicle_global_id getVes2idx() const;

    linkspec_global_id getLinkSpec1gidx() const;
    linkspec_global_id getLinkSpec2gidx() const;

    // The products
    spec_global_id getSpec1gidx() const;
    spec_global_id getSpec2gidx() const;

    // Placeholder. Might need to support other orders at some stage
    inline uint order() const noexcept {
        return 1;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    vesunbind_global_id pIdx;
    std::string pName;
    double pKcst;

    int pImmobility;

    // Need to copy this data for setup
    std::pair<model::Vesicle*, model::LinkSpec*> pLinks1;
    std::pair<model::Vesicle*, model::LinkSpec*> pLinks2;
    std::pair<model::Vesicle*, model::Spec*> pProducts1;
    std::pair<model::Vesicle*, model::Spec*> pProducts2;

    vesicle_global_id pVesicle_1_idx;
    vesicle_global_id pVesicle_2_idx;

    linkspec_global_id pLinkSpec_1_gidx;
    linkspec_global_id pLinkSpec_2_gidx;

    spec_global_id pSpec_1_gidx;
    spec_global_id pSpec_2_gidx;

    bool pSetupdone;
};

}  // namespace steps::solver
