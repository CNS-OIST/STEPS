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
#include "model/raftgen.hpp"
#include "model/spec.hpp"
#include "solver/api.hpp"
#include "solver/fwd.hpp"
#include "solver/raftdef.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

/// Defined Raft Genesis
class RaftGendef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the raft genesis reaction.
    /// \param endo Pointer to the RaftGen object.
    RaftGendef(Statedef* sd, raftgen_global_id idx, model::RaftGen* raftgen);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT GENESIS RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this raft genesis rule.
    inline raftgen_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the raft genesis.
    inline std::string const& name() const noexcept {
        return pName;
    }

    inline double kcst() const noexcept {
        return pKcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    //////////////////////// Membrane species //////////////////////////////

    /// Returns a description of how an occurence of this raft genesis
    /// depends on some species, defined by its global index idx, to occur.
    depT dep_S(spec_global_id gidx) const;

    bool reqspec_S(spec_global_id gidx) const;

    inline uint countSpecs_S() const noexcept {
        return statedef()->countSpecs();
    }

    const auto& lhs_S() const noexcept {
        return pSpec_S_LHS;
    }

    // For clarity
    inline uint countSpecs_global() const noexcept {
        return statedef()->countSpecs();
    }

    ////////////////////////////////////////////////////////////////////////

    inline Statedef* statedef() const noexcept {
        return pStatedef;
    }

    inline Raftdef* raftdef() const noexcept {
        return pRaftdef;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    raftgen_global_id pIdx;

    std::string pName;

    double pKcst;

    bool pSetupdone;

    // Only used during setup
    model::SpecPVec pSDeps;
    model::Raft* pRaft;

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

    // Store the raft of creation as a Raftdef pointer
    Raftdef* pRaftdef{};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver