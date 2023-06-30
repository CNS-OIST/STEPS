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
#include "model/linkspec.hpp"
#include "model/spec.hpp"
#include "model/vesbind.hpp"
#include "model/vesicle.hpp"
#include "solver/api.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

/// Defined vesicle-vesicle binding reaction
class VesBinddef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the vesicle binding reaction.
    /// \param vb Pointer to the associated VesBind object.
    VesBinddef(Statedef* sd, vesbind_global_id idx, model::VesBind* vb);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESBIND RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this vesicle binding reaction rule.
    inline vesbind_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the vesicle binding reaction.
    const std::string& name() const;

    /// Return the MACROscopic vesicle binding reaction constant.
    double kcst() const;

    inline double max_distance() const noexcept {
        return pMaxDistance;
    }

    inline double min_distance() const noexcept {
        return pMinDistance;
    }

    // Added for convenience for solver object access
    inline uint countSpecsGlobal() const noexcept {
        return pStatedef->countSpecs();
    }

    inline uint countLinkSpecsGlobal() const noexcept {
        return pStatedef->countLinkSpecs();
    }

    int immobility() const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    vesicle_global_id getVes1idx() const;
    vesicle_global_id getVes2idx() const;

    spec_global_id getSpec1gidx() const;
    spec_global_id getSpec2gidx() const;

    linkspec_global_id getLinkSpec1gidx() const;
    linkspec_global_id getLinkSpec2gidx() const;

    LinkSpecdef* getLinkSpec1def() const;
    LinkSpecdef* getLinkSpec2def() const;

    uint vdep1(spec_global_id gidx) const;
    uint vdep2(spec_global_id gidx) const;

    uint ldep1(linkspec_global_id gidx) const;
    uint ldep2(linkspec_global_id gidx) const;

    // Placeholder. Might need to support other orders at some stage
    inline uint order() const noexcept {
        return 2;
    }

    bool gotVDep1() const noexcept {
        return pGotVDep1;
    }
    bool gotVDep2() const noexcept {
        return pGotVDep2;
    }

    bool gotLDep1() const noexcept {
        return pGotLDep1;
    }
    bool gotLDep2() const noexcept {
        return pGotLDep2;
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
    vesbind_global_id pIdx;
    std::string pName;
    double pKcst;

    double pMaxDistance;
    double pMinDistance;

    // Need to copy this data for setup
    std::pair<model::Vesicle*, model::Spec*> pReactants1;
    std::pair<model::Vesicle*, model::Spec*> pReactants2;

    std::pair<model::Vesicle*, model::LinkSpec*> pProducts1;
    std::pair<model::Vesicle*, model::LinkSpec*> pProducts2;

    LinkSpecdef* pProduct1def{};
    LinkSpecdef* pProduct2def{};

    vesicle_global_id pVesicle_1_idx;
    vesicle_global_id pVesicle_2_idx;
    spec_global_id pSpec_1_gidx;
    spec_global_id pSpec_2_gidx;
    linkspec_global_id pProduct1_gidx;
    linkspec_global_id pProduct2_gidx;

    model::SpecPVec pVdep1;
    model::SpecPVec pVdep2;
    model::LinkSpecPVec pLdep1;
    model::LinkSpecPVec pLdep2;

    int pImmobility;

    bool pSetupdone;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    // Information about the vesicle surface species dependency
    util::strongid_vector<spec_global_id, uint> pSpec_VDEP1;
    util::strongid_vector<spec_global_id, uint> pSpec_VDEP2;
    util::strongid_vector<linkspec_global_id, uint> pSpec_LDEP1;
    util::strongid_vector<linkspec_global_id, uint> pSpec_LDEP2;

    bool pGotVDep1;
    bool pGotVDep2;
    bool pGotLDep1;
    bool pGotLDep2;
};

}  // namespace steps::solver
