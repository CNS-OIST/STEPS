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

#include <gsl>

// STEPS headers.
#include "model/raft.hpp"
#include "solver/fwd.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

// logging
#include "easylogging++.h"

namespace steps::solver {

/// Defined Raft
class Raftdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the raft.
    /// \param r Pointer to the associated Raft object.
    Raftdef(Statedef* sd, raft_global_id idx, model::Raft* r);

    /// Destructor
    ~Raftdef();

    void reset() const;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this raft.
    inline raft_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the raft.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the diameter of the raft.
    double diameter() const noexcept {
        return pDiameter;
    }

    /// Return the diffusion constant.
    inline double dcst() const noexcept {
        return pDcst;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    void setup_references();
    void setup_indices();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of global species
    // Note: helper function for tetvesicle::Raft so not private
    inline uint countSpecs_global() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for this raft.
    // Note: used by tetvesicle::Raft
    inline uint countSpecs() const noexcept {
        return pSpecsN_Rs;
    }

  private:
    /// Return the number of species defined for the outer compartment.
    // For rafts we deal with global indices since many different types of
    // compartments could be used
    inline uint countSpecs_O() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for the surface (patch).
    // For rafts we deal with global indices since many different types of
    // patches could be used
    inline uint countSpecs_S() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for the inner compartment.
    // For rafts we deal with global indices since many different types of
    // compartments could be used
    inline uint countSpecs_I() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for the raft surface dependency.
    inline uint countSpecs_RsDep() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for the raft surface anti-dependency.
    inline uint countSpecs_AntiRsDep() const noexcept {
        return pSpecsN_global;
    }

  public:
    /// Return the local species index for global index argument.
    ///
    /// \param gidx Global index of the species.
    inline spec_local_id specG2L(spec_global_id gidx) const noexcept {
        return pSpec_G2L[gidx.get()];
    }

    /// Return the global species index for local index argument.
    ///
    /// \param lidx local index of the species.
    inline spec_global_id specL2G(spec_local_id lidx) const noexcept {
        return pSpec_L2G[lidx.get()];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft surface reactions that can occur on this
    /// raft.
    inline uint countRaftSReacs() const noexcept {
        return pRaftSReacsN;
    }

    /// Return the local raft surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    inline raftsreac_local_id raftsreacG2L(raftsreac_global_id gidx) const noexcept {
        return pRaftSReac_G2L[gidx.get()];
    }

    /// Return the global raft surface reaction index for global index argument.
    ///
    /// \param gidx Local index of the raft surface reaction.
    inline raftsreac_global_id raftsreacL2G(raftsreac_local_id lidx) const noexcept {
        return pRaftSReac_L2G[lidx.get()];
    }

    /// Return a pointer to reaction definition object (type RaftSReacdef)
    /// specified by local index.
    RaftSReacdef* raftsreacdef(raftsreac_local_id lidx) const;

    int raftsreac_dep_Rs(raftsreac_local_id srlidx, spec_local_id splidx) const noexcept;
    // O uses global species arguments
    int raftsreac_dep_O(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept;
    // S uses global species arguments
    int raftsreac_dep_S(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept;
    // I uses global species arguments
    int raftsreac_dep_I(raftsreac_local_id srlidx, spec_global_id spgidx) const noexcept;

  private:
    /// internal utility function

    template <typename T>
    static constexpr gsl::span<const T> to_span(const std::vector<T>& container,
                                                raftsreac_local_id lidx,
                                                uint size) {
        const auto begin = container.data() + lidx.get() * size;
        return {begin, begin + size};
    }

    template <typename StrongId, typename T>
    static constexpr auto to_strong_span(const std::vector<T>& container,
                                         raftsreac_local_id lidx,
                                         uint size) {
        return util::make_strong_random_accessor<StrongId>(to_span(container, lidx, size));
    }


  public:
    ////////////////////////////////////////////////////////////////////////
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.

    strongid_span<spec_local_id, const uint> raftsreac_lhs_Rs(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pRaftSReac_LHS_Rs_Spec, lidx, countSpecs());
    }

    strongid_span<spec_global_id, const uint> raftsreac_lhs_O(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_LHS_O_Spec, lidx, countSpecs_O());
    }

    strongid_span<spec_global_id, const uint> raftsreac_lhs_S(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_LHS_S_Spec, lidx, countSpecs_S());
    }

    strongid_span<spec_global_id, const uint> raftsreac_lhs_I(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_LHS_I_Spec, lidx, countSpecs_I());
    }

    // Raft surface dependencies
    strongid_span<spec_global_id, const uint> raftsreac_rsdep(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_RsDEP_Spec, lidx, countSpecs_RsDep());
    }

    strongid_span<spec_global_id, const uint> raftsreac_anti_rsdep(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_AntiRsDEP_Spec,
                                              lidx,
                                              countSpecs_AntiRsDep());
    }

    ////////////////////////////////////////////////////////////////////////
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.

    strongid_span<spec_local_id, const int> raftsreac_upd_Rs(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pRaftSReac_UPD_Rs_Spec, lidx, countSpecs());
    }

    strongid_span<spec_global_id, const int> raftsreac_upd_O(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_UPD_O_Spec, lidx, countSpecs_O());
    }

    strongid_span<spec_global_id, const int> raftsreac_upd_S(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_UPD_S_Spec, lidx, countSpecs_S());
    }

    strongid_span<spec_global_id, const int> raftsreac_upd_I(
        raftsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pRaftSReac_UPD_I_Spec, lidx, countSpecs_I());
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: ENDOCYTOTIC REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of endocytotic reactions that can occur on this
    /// raft.
    inline uint countRaftEndocytosis() const noexcept {
        return pRaftEndocytosisN;
    }

    /// Return the local endocytotic reaction index for global index argument.
    ///
    /// \param gidx Global index of the endocytotic reaction.
    inline raftendocytosis_local_id raftendocytosisG2L(
        raftendocytosis_global_id gidx) const noexcept {
        return pRaftEndocytosis_G2L[gidx.get()];
    }

    /// Return the global endocytotic reaction index for local index argument.
    ///
    /// \param lidx Local index of the endocytotic reaction.
    inline raftendocytosis_global_id raftendocytosisL2G(
        raftendocytosis_local_id lidx) const noexcept {
        return pRaftEndocytosis_L2G[lidx.get()];
    }

    /// Return a pointer to reaction definition object (type RaftEndocytosis)
    /// specified by local index.
    RaftEndocytosisdef* raftendocytosisdef(raftendocytosis_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT DISSOLUTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft dissolution reactions that can occur on
    /// this patch.
    inline uint countRaftDiss() const noexcept {
        return pRaftDisN;
    }

    /// Return the local raft dissolution index for global index argument.
    ///
    /// \param gidx Global index of the raft dissolution
    inline raftdis_local_id raftdisG2L(raftdis_global_id gidx) const noexcept {
        return pRaftDis_G2L[gidx.get()];
    }

    /// Return the global raft dissolution index for local index argument.
    ///
    /// \param gidx Local index of the raft dissolution
    inline raftdis_global_id raftdisL2G(raftdis_local_id lidx) const noexcept {
        return pRaftDis_L2G[lidx.get()];
    }

    /// Return a pointer to reaction definition object (type RaftDisdef)
    /// specified by local index.
    RaftDisdef* raftdisdef(raftdis_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    Statedef* pStatedef;
    raft_global_id pIdx;
    std::string pName;
    double pDiameter;

    // The diffusion constant
    double pDcst;

    // Keep track of whether setup methods have been called.
    bool pSetupRefsdone;
    bool pSetupIndsdone;

    // The enclosed surface systems, stored as strings
    std::set<std::string> pRssys;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Number of species embedded in raft surface (_Rs)
    ///
    uint pSpecsN_Rs;
    // Note, pSpecsN_S, _I and _O are not needed because we use global indices for
    // those

    // Multiple 'outer' volumes can be supported so just store the
    // global number of species
    uint pSpecsN_global;

    /// Table to resolve species index (global -> local).
    std::vector<spec_local_id> pSpec_G2L;

    /// Table to resolve species index (local -> global).
    std::vector<spec_global_id> pSpec_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: RAFT SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of surface reactions occurring in patch.
    uint pRaftSReacsN;

    /// Table to resolve reaction rule indices (global -> local).
    std::vector<raftsreac_local_id> pRaftSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    std::vector<raftsreac_global_id> pRaftSReac_L2G;

    // Species indices are local only on the raft surface
    inline uint _IDX_RaftSReac_Rs_Spec(raftsreac_local_id srlidx) const noexcept {
        return countSpecs() * srlidx.get();
    }
    inline uint _IDX_RaftSReac_Rs_Spec(raftsreac_local_id srlidx,
                                       spec_local_id splidx) const noexcept {
        return (countSpecs() * srlidx.get()) + splidx.get();
    }
    inline uint _IDX_RaftSReac_O_Spec_global(raftsreac_local_id srlidx) const noexcept {
        return countSpecs_O() * srlidx.get();
    }
    inline uint _IDX_RaftSReac_O_Spec_global(raftsreac_local_id srlidx,
                                             spec_global_id spgidx) const noexcept {
        return (countSpecs_O() * srlidx.get()) + spgidx.get();
    }
    inline uint _IDX_RaftSReac_S_Spec_global(raftsreac_local_id srlidx) const noexcept {
        return countSpecs_S() * srlidx.get();
    }
    inline uint _IDX_RaftSReac_S_Spec_global(raftsreac_local_id srlidx,
                                             spec_global_id spgidx) const noexcept {
        return (countSpecs_S() * srlidx.get()) + spgidx.get();
    }
    inline uint _IDX_RaftSReac_I_Spec_global(raftsreac_local_id srlidx) const noexcept {
        return countSpecs_I() * srlidx.get();
    }
    inline uint _IDX_RaftSReac_I_Spec_global(raftsreac_local_id srlidx,
                                             spec_global_id spgidx) const noexcept {
        return (countSpecs_I() * srlidx.get()) + spgidx.get();
    }

    inline uint _IDX_RaftSReac_Rsdep_Spec_global(raftsreac_local_id srlidx,
                                                 spec_global_id spgidx) const noexcept {
        return (countSpecs_RsDep() * srlidx.get()) + spgidx.get();
    }

    inline uint _IDX_RaftSReac_AntiRsdep_Spec_global(raftsreac_local_id srlidx,
                                                     spec_global_id spgidx) const noexcept {
        return (countSpecs_AntiRsDep() * srlidx.get()) + spgidx.get();
    }

    std::vector<int> pRaftSReac_DEP_Rs_Spec;
    std::vector<int> pRaftSReac_DEP_O_Spec;
    std::vector<int> pRaftSReac_DEP_S_Spec;
    std::vector<int> pRaftSReac_DEP_I_Spec;
    std::vector<uint> pRaftSReac_LHS_Rs_Spec;
    std::vector<uint> pRaftSReac_LHS_O_Spec;
    std::vector<uint> pRaftSReac_LHS_S_Spec;
    std::vector<uint> pRaftSReac_LHS_I_Spec;
    std::vector<int> pRaftSReac_UPD_Rs_Spec;
    std::vector<int> pRaftSReac_UPD_O_Spec;
    std::vector<int> pRaftSReac_UPD_S_Spec;
    std::vector<int> pRaftSReac_UPD_I_Spec;

    std::vector<uint> pRaftSReac_RsDEP_Spec;
    std::vector<uint> pRaftSReac_AntiRsDEP_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: RAFT SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of endocytosis occurring in raft.
    uint pRaftEndocytosisN;

    /// Table to resolve endocytosis rule indices (global -> local).
    std::vector<raftendocytosis_local_id> pRaftEndocytosis_G2L;

    /// Table to resolve endocytosis rule indices (local -> global).
    std::vector<raftendocytosis_global_id> pRaftEndocytosis_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: RAFT DISSOLUTION
    ////////////////////////////////////////////////////////////////////////

    /// Number of raft dissolution occurring in raft.
    uint pRaftDisN;

    /// Table to resolve raft dissolution rule indices (global -> local).
    std::vector<raftdis_local_id> pRaftDis_G2L;

    /// Table to resolve raft dissolution rule indices (local -> global).
    std::vector<raftdis_global_id> pRaftDis_L2G;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::solver
