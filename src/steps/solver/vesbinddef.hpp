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

#include <iosfwd>
#include <string>

#include "fwd.hpp"
#include "model/fwd.hpp"
#include "model/vesbind.hpp"

namespace steps::solver {

/// Defined vesicle-vesicle binding reaction
class VesBinddef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the vesicle binding reaction.
    /// \param vb Reference to the associated VesBind object.
    VesBinddef(Statedef& sd, vesbind_global_id idx, model::VesBind& vb);

    VesBinddef(const VesBinddef&) = delete;
    VesBinddef& operator=(const VesBinddef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

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
    inline const std::string& name() const noexcept {
        return pName;
    }

    /// Return the MACROscopic vesicle binding reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    inline double max_distance() const noexcept {
        return pMaxDistance;
    }

    inline double min_distance() const noexcept {
        return pMinDistance;
    }

    // Added for convenience for solver object access
    inline uint countSpecsGlobal() const noexcept {
        return pSpec_VDEP1.size();
    }

    inline uint countLinkSpecsGlobal() const noexcept {
        return pSpec_LDEP1.size();
    }

    inline int immobility() const noexcept {
        return pImmobility;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    vesicle_global_id getVes1idx() const;
    vesicle_global_id getVes2idx() const;

    spec_global_id getSpec1gidx() const;
    spec_global_id getSpec2gidx() const;

    linkspec_global_id getLinkSpec1gidx() const;
    linkspec_global_id getLinkSpec2gidx() const;

    LinkSpecdef& getLinkSpec1def() const;
    LinkSpecdef& getLinkSpec2def() const;

    uint vdep1(spec_global_id gidx) const;
    uint vdep2(spec_global_id gidx) const;

    uint ldep1(linkspec_global_id gidx) const;
    uint ldep2(linkspec_global_id gidx) const;

    // Placeholder. Might need to support other orders at some stage
    inline uint order() const noexcept {
        return 2;
    }

    bool gotVDep1() const noexcept {
        return !pVdep1.empty();
    }
    bool gotVDep2() const noexcept {
        return !pVdep2.empty();
    }

    bool gotLDep1() const noexcept {
        return !pLdep1.empty();
    }
    bool gotLDep2() const noexcept {
        return !pLdep2.empty();
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    void setup(const Statedef& sd);

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    const vesbind_global_id pIdx;
    const std::string pName;
    const double pKcst;

    const double pMaxDistance;
    const double pMinDistance;

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

    const std::vector<model::Spec*> pVdep1;
    const std::vector<model::Spec*> pVdep2;
    const std::vector<model::LinkSpec*> pLdep1;
    const std::vector<model::LinkSpec*> pLdep2;

    const int pImmobility;

    bool pSetupdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    // Information about the vesicle surface species dependency
    util::strongid_vector<spec_global_id, uint> pSpec_VDEP1;
    util::strongid_vector<spec_global_id, uint> pSpec_VDEP2;
    util::strongid_vector<linkspec_global_id, uint> pSpec_LDEP1;
    util::strongid_vector<linkspec_global_id, uint> pSpec_LDEP2;
};

}  // namespace steps::solver
