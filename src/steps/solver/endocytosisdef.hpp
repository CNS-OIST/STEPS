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
#include <vector>

#include "fwd.hpp"
#include "model/endocytosis.hpp"
#include "model/fwd.hpp"
#include "statedef.hpp"
#include "util/vocabulary.hpp"

namespace steps::solver {

// Events struct
class EndocytosisEvent {
  public:
    EndocytosisEvent()
        : time(0)
        , tidx(triangle_global_id::unknown_value())
        , vidx(vesicle_individual_id::unknown_value()) {}
    EndocytosisEvent(double t, triangle_global_id ti, vesicle_individual_id v)
        : time(t)
        , tidx(ti)
        , vidx(v) {}
    double time;
    triangle_global_id tidx;
    vesicle_individual_id vidx;
};

/// Defined Endocytisis Reaction.
class Endocytosisdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the endocytotic reaction.
    /// \param endo Reference to the Endocytosis object.
    Endocytosisdef(Statedef& sd, endocytosis_global_id idx, model::Endocytosis& endo);

    Endocytosisdef(const Endocytosisdef&) = delete;
    Endocytosisdef& operator=(const Endocytosisdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: ENDOCYTOTIC REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this endocytotic reaction rule.
    inline endocytosis_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the endocytotic reaction.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the MACROscopic reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup(const Statedef&);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    //////////////////////// Membrane species //////////////////////////////

    /// Returns a description of how an occurence of this endocytosis
    /// depends on some species, defined by its global index idx, to occur.
    depT dep_S(spec_global_id gidx) const;

    bool reqspec_S(spec_global_id gidx) const;

    inline uint countSpecs_S() const noexcept {
        return pStatedef.countSpecs();
    }

    inline const auto& lhs_S() const noexcept {
        return pSpec_S_LHS;
    }

    //////////////////////////// Vesicles //////////////////////////////////

    /// Returns the defined vesicle product.
    ///
    Vesicledef& rhs_I_ves() const;

    /// Returns the idx of the vesicle product.
    ///
    vesicle_global_id rhs_I_ves_uint() const;

    inline bool inner() const noexcept {
        return pInner;
    }

    ////////////////////////////////////////////////////////////////////////

    // May be prudent to store pointer to vesicledef, but for now let's
    // do it this longwinded way

    inline const Statedef& statedef() const noexcept {
        return pStatedef;
    }

  private:
    Statedef& pStatedef;
    const endocytosis_global_id pIdx;

    const std::string pName;

    const double pKcst;

    // The stoichiometry stored as model level Spec objects.
    // To be used during setup ONLY
    const model::Vesicle& pIrhs;

    bool pSetupdone{false};

    // Only used during setup
    std::vector<model::Spec*> pSDeps;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    // The following two are used like usual STEPS arrays to say how this
    // SSA event depends on the presence of certain species
    util::strongid_vector<spec_global_id, uint> pSpec_S_LHS;
    util::strongid_vector<spec_global_id, depT> pSpec_S_DEP;

    Vesicledef* pVes_I_RHS{nullptr};
    vesicle_global_id pVes_I_RHS_id;

    bool pInner;
};

}  // namespace steps::solver
