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
#include "model/raftendocytosis.hpp"
#include "solver/fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps::solver {

// Events struct
class RaftEndocytosisEvent {
  public:
    RaftEndocytosisEvent()
        : time(0)
        , ridx(raft_individual_id::unknown_value())
        , tidx(triangle_global_id::unknown_value())
        , vidx(vesicle_individual_id::unknown_value()) {}
    RaftEndocytosisEvent(double t,
                         raft_individual_id r,
                         triangle_global_id ti,
                         vesicle_individual_id v)
        : time(t)
        , ridx(r)
        , tidx(ti)
        , vidx(v) {}
    double time;
    raft_individual_id ridx;
    triangle_global_id tidx;
    vesicle_individual_id vidx;
};

/// Defined RaftEndocytisis Reaction.
class RaftEndocytosisdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the raft endocytotic reaction.
    /// \param endo Reference to the RaftEndocytosis object.
    RaftEndocytosisdef(Statedef& sd, raftendocytosis_global_id idx, model::RaftEndocytosis& endo);

    RaftEndocytosisdef(const RaftEndocytosisdef&) = delete;
    RaftEndocytosisdef& operator=(const RaftEndocytosisdef&) = delete;

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
    inline raftendocytosis_global_id gidx() const noexcept {
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
    // SOLVER METHODS
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup(const Statedef& sd);

    void reset();

    void setKcst(double kcst);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    //////////////////////// Membrane species //////////////////////////////

    /// Returns a description of how an occurence of this raft endocytosis
    /// reaction depends on some species, defined by its global index idx, to
    /// occur.
    depT dep_S(spec_global_id gidx) const;

    bool reqspec_S(spec_global_id gidx) const;

    inline uint countSpecs_S() const noexcept {
        return pSpec_S_DEP.size();
    }
    inline const auto& lhs_S() const noexcept {
        return pSpec_S_LHS;
    }

    // Just for clarity (instead of using countSpecs_S)- surface species count
    // might change in the future
    inline uint countSpecs_global() const noexcept {
        return pCountSpecs;
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

    ///////////////////////////// Extent ////////////////////////////////////

    inline unsigned long getExtent() const noexcept {
        return pExtent;
    }

    std::vector<RaftEndocytosisEvent> getEvents();

    inline void addEvent(double time,
                         raft_individual_id ridx,
                         triangle_global_id tidx,
                         vesicle_individual_id vidx) noexcept {
        pExtent++;
        pEvents.emplace_back(time, ridx, tidx, vidx);
    }
    ////////////////////////////////////////////////////////////////////////

    // May be prudent to store pointer to vesicledef, but for now let's
    // do it this longwinded way

    inline const Statedef& statedef() const noexcept {
        return pStatedef;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    const Statedef& pStatedef;
    const raftendocytosis_global_id pIdx;

    const std::string pName;

    double pKcst;

    // Store now the default Kcst because we're going to allow a global change
    // by an API call setRaftEndocytosisK
    const double pDefaultKcst;

    // The stoichiometry stored as model level Spec objects.
    // To be used during setup ONLY
    const model::Vesicle& pIrhs;

    const uint pCountSpecs;

    bool pSetupdone{false};

    // Only used during setup
    const std::vector<model::Spec*> pSDeps;

    unsigned long pExtent;

    std::vector<RaftEndocytosisEvent> pEvents;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    // The following two are used like usual STEPS arrays to say how this
    // SSA event depends on the presence of certain species
    util::strongid_vector<spec_global_id, uint> pSpec_S_LHS;
    // Dep is going to look a bit different here from dep in the model object- LHS
    // will store possible multiple copies of a species (like 'dep' in the model),
    // DEP here will be a yes or no
    util::strongid_vector<spec_global_id, depT> pSpec_S_DEP;

    // We've stopped using strings- let's use vesdef pointer instead
    Vesicledef* pVes_I_RHS;

    vesicle_global_id pVes_I_RHS_uint;

    const bool pInner;
};

}  // namespace steps::solver
