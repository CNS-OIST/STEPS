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
#include "model/exocytosis.hpp"
#include "model/raft.hpp"
#include "model/spec.hpp"
#include "model/vesicle.hpp"
#include "solver/api.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

// Events struct
class ExocytosisEvent {
  public:
    ExocytosisEvent()
        : time(0)
        , vidx(vesicle_individual_id::unknown_value())
        , tidx(triangle_global_id::unknown_value())
        , ridx(raft_individual_id::unknown_value()) {}
    ExocytosisEvent(double t, vesicle_individual_id v, triangle_global_id ti, raft_individual_id r)
        : time(t)
        , vidx(v)
        , tidx(ti)
        , ridx(r) {}
    double time;
    vesicle_individual_id vidx;
    triangle_global_id tidx;
    raft_individual_id ridx;
};

/// Defined Exocytosis Reaction.
class Exocytosisdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the exocytotic reaction.
    /// \param endo Pointer to the Exocytosis object.
    Exocytosisdef(Statedef* sd, exocytosis_global_id idx, model::Exocytosis* exo);

    /// Destructor
    ~Exocytosisdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EXOCYTOTIC REACTION RULE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this exocytotic reaction rule.
    inline exocytosis_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the exocytotic reaction.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the MACROscopic reaction constant.
    inline double kcst() const noexcept {
        return pKcst;
    }

    void setKcst(double k);

    // Need reset function to reset extent
    inline void reset() {
        pKcst = pExocytosis->getKcst();
        pExtent = 0;
        pEvents.clear();
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    /// Returns the number of molecules of species idx required in
    /// the surface patch (_S)
    /// to activate exocytosis.
    ///
    uint lhs_V(spec_global_id gidx) const;

    const auto& lhs_V() const noexcept {
        return pSpec_V_LHS;
    }

    uint countSpecs_V() const noexcept {
        return statedef()->countSpecs();
    }

    uint vdepssize() const noexcept {
        return pVDeps.size();
    }

    ////////////////////////////////////////////////////////////////////////

    /// Returns a description of how an occurrence of this exocytosis reaction
    /// depends on some surface species, defined by its global index idx, to
    /// occur.
    ///
    depT dep_V(spec_global_id gidx) const;

    inline Statedef* statedef() const noexcept {
        return pStatedef;
    }

    inline Raftdef* raftdef() const noexcept {
        return pRaftdef;
    }

    inline bool getKissAndRun() const noexcept {
        return pKissAndRun;
    }

    inline const std::map<spec_global_id, spec_global_id>& getKissAndRunSpecChanges()
        const noexcept {
        return pKissAndRunSpecChanges;
    }

    inline unsigned long getExtent() const noexcept {
        return pExtent;
    }

    std::vector<ExocytosisEvent> getEvents();

    inline void addEvent(double time,
                         vesicle_individual_id vidx,
                         triangle_global_id tidx,
                         raft_individual_id ridx) noexcept {
        pExtent++;
        pEvents.emplace_back(time, vidx, tidx, ridx);
    }

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    exocytosis_global_id pIdx;

    std::string pName;

    double pKcst;

    // Have to store this to access the original kcsts to allow for reset
    model::Exocytosis* pExocytosis;

    bool pSetupdone;

    unsigned long pExtent;

    std::vector<ExocytosisEvent> pEvents;

    ////////////////////////////////////////////////////////////////////////
    // DATA: STOICHIOMETRY
    ////////////////////////////////////////////////////////////////////////

    // Only used during setup
    model::SpecPVec pVDeps;
    model::Raft* pRaft;

    // The following two are used like usual STEPS arrays to say how this
    // SSA event depends on the presence of certain species
    util::strongid_vector<spec_global_id, uint> pSpec_V_LHS;
    util::strongid_vector<spec_global_id, depT> pSpec_V_DEP;

    // Store the raft of creation as a Raftdef pointer
    Raftdef* pRaftdef;

    // Store whether this is a kiss-n-run exocytosis event
    bool pKissAndRun;

    // Store any species that will be transfered from vesicle to membrane upon kiss-n-run
    std::map<spec_global_id, spec_global_id> pKissAndRunSpecChanges;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver