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

#include "model/fwd.hpp"
#include "model/vesicle.hpp"
#include "solver/fwd.hpp"

#include <gsl>

namespace steps::solver {

/// Defined Vesicle
class Vesicledef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the vesicle.
    /// \param d Reference to the associated Vesicle object.
    Vesicledef(Statedef& sd, vesicle_global_id idx, model::Vesicle& d);

    Vesicledef(const Vesicledef&) = delete;
    Vesicledef& operator=(const Vesicledef&) = delete;

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this vesicle.
    inline vesicle_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the vesicle.
    inline std::string const& name() const noexcept {
        return pName;
    }

    /// Return the diameter of the vesicle.
    inline double diameter() const noexcept {
        return pDiameter;
    }

    /// Return the diffusion constant.
    inline double dcst() const noexcept {
        return pDcst;
    }

    // Return the inner volume of the vesicle
    double vol() const noexcept;

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

  private:
    /// Return the number of species defined for this vesicle.
    inline uint countSpecs() const noexcept {
        return pSpecsN_V;
    }

    /* Not needed because link specs use global indices??
    /// Return the number of species defined for this vesicle.
    inline uint countLinkSpecs(void) const
    { return pSpecsN_L; }
    */

    /// Return the number of species defined for the outer compartment.
    // For vesicles we deal with global indices since many different types of
    // compartments could be used
    inline uint countSpecs_O() const noexcept {
        return pSpecsN_global;
    }

    /// Return the number of species defined for the surface (patch).
    // For vesicles we deal with global indices since many different types of
    // patches could be used
    inline uint countSpecs_S() const noexcept {
        return pSpecsN_global;
    }

    inline uint countSpecs_VDep() const noexcept {
        return pSpecsN_global;
    }

  public:
    // Note: countSpecs_I() is used by vesproxy
    /// Return the number of species defined for the inner compartment.
    // For vesicles we deal with global indices since many different types of
    // compartments could be used
    inline uint countSpecs_I() const noexcept {
        return pSpecsN_global;
    }

    /// Return the local species index for global index argument.
    ///
    /// \param gidx Global index of the species.
    inline spec_local_id specG2L(spec_global_id gidx) const noexcept {
        return pSpec_G2L[gidx];
    }

    /// Return the global species index for local index argument.
    ///
    /// \param lidx local index of the species.
    inline spec_global_id specL2G(spec_local_id lidx) const noexcept {
        return pSpec_L2G[lidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: LINK SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of link species defined for this vesicle.
    // We'll just use global indices because that's how they are stored in
    // LinkSpec objects
    inline uint countLinkSpecs_V() const noexcept {
        return pLinkSpecsN_global;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface reactions that can occur on this
    /// vesicle.
    inline uint countVesSReacs() const noexcept {
        return pVesSReacsN;
    }

    /// Return the local surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    inline vessreac_local_id vessreacG2L(vessreac_global_id gidx) const noexcept {
        return pVesSReac_G2L[gidx];
    }

    /// Return the global surface reaction index for local index argument.
    ///
    /// \param gidx Local index of the surface reaction.
    inline vessreac_global_id vessreacL2G(vessreac_local_id lidx) const noexcept {
        return pVesSReac_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type VesSReacdef)
    /// specified by local index.
    VesSReacdef& vessreacdef(vessreac_local_id lidx) const;

    // depV uses local species arguments
    int vessreac_dep_V(vessreac_local_id, spec_local_id) const noexcept;

    // L uses global arguments
    int vessreac_dep_L(vessreac_local_id, linkspec_global_id) const noexcept;
    // O uses global species arguments
    int vessreac_dep_O(vessreac_local_id, spec_global_id) const noexcept;
    // S uses global species arguments
    int vessreac_dep_S(vessreac_local_id, spec_global_id) const noexcept;

  private:
    /// internal utility function

    template <typename T>
    static constexpr gsl::span<const T> to_span(const std::vector<T>& container,
                                                vessreac_local_id lidx,
                                                uint size) {
        const auto begin = container.data() + lidx.get() * size;
        return {begin, begin + size};
    }

    template <typename StrongId, typename T>
    static constexpr auto to_strong_span(const std::vector<T>& container,
                                         vessreac_local_id lidx,
                                         uint size) {
        return util::make_strong_random_accessor<StrongId>(to_span(container, lidx, size));
    }

  public:
    ///
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const uint> vessreac_lhs_V(vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVesSReac_LHS_V_Spec, lidx, countSpecs());
    }

    strongid_span<linkspec_global_id, const uint> vessreac_lhs_L(
        vessreac_local_id lidx) const noexcept {
        return to_strong_span<linkspec_global_id>(pVesSReac_LHS_L_Spec, lidx, countLinkSpecs_V());
    }

    strongid_span<spec_global_id, const uint> vessreac_lhs_O(
        vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_LHS_O_Spec, lidx, countSpecs_O());
    }

    strongid_span<spec_global_id, const uint> vessreac_lhs_S(
        vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_LHS_S_Spec, lidx, countSpecs_S());
    }

    // Vesicle surface dependencies
    strongid_span<spec_global_id, const uint> vessreac_vdep(vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_VDEP_Spec, lidx, countSpecs_VDep());
    }

    ///
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const int> vessreac_upd_V(vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVesSReac_UPD_V_Spec, lidx, countSpecs());
    }

    // L uses global species arguments
    strongid_span<linkspec_global_id, const int> vessreac_upd_L(
        vessreac_local_id lidx) const noexcept {
        return to_strong_span<linkspec_global_id>(pVesSReac_UPD_L_Spec, lidx, countLinkSpecs_V());
    }

    // O uses global species arguments
    strongid_span<spec_global_id, const int> vessreac_upd_O(vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_UPD_O_Spec, lidx, countSpecs_O());
    }

    // S uses global species arguments
    strongid_span<spec_global_id, const int> vessreac_upd_S(vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_UPD_S_Spec, lidx, countSpecs_S());
    }

    // I uses global species arguments. Special case that update is
    // always positive because it can't be consumed
    strongid_span<spec_global_id, const uint> vessreac_upd_I(
        vessreac_local_id lidx) const noexcept {
        return to_strong_span<spec_global_id>(pVesSReac_UPD_I_Spec, lidx, countSpecs_I());
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EXOCYTOTIC REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    static const uint INACTIVATED = 1;

    /// Return the total number of exocytotic reactions that can occur on this
    /// patch.
    inline uint countExocytosis() const noexcept {
        return pExocytosisN;
    }

    /// Return the local exocytotic reaction index for global index argument.
    ///
    /// \param gidx Global index of the exocytotic reaction.
    inline exocytosis_local_id exocytosisG2L(exocytosis_global_id gidx) const noexcept {
        return pExocytosis_G2L[gidx];
    }

    /// Return the global exocytotic reaction index for local index argument.
    ///
    /// \param lidx Local index of the exocytotic reaction.
    inline exocytosis_global_id exocytosisL2G(exocytosis_local_id lidx) const noexcept {
        return pExocytosis_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type Exocytosisdef)
    /// specified by local index.
    Exocytosisdef& exocytosisdef(exocytosis_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: EXOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the kcst for exocytotic reaction specified by local index
    ///
    /// \param lidx Local index of the exocytotic reaction.
    /// \param kcst Rate constant of the exocytotic reaction.
    void setExoKcst(exocytosis_local_id lidx, double kcst);

    /// Activate or inactivate a exocytotic reaction specified by local index.
    ///
    /// \param lidx Local index of the exocytotic reaction.
    /// \param active Flag to activate / inactivate the exocytotic reaction.
    void setExoActive(exocytosis_local_id lidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules for this vesicle.
    inline uint countVesSurfDiffs() const noexcept {
        return pVesSDiffsN;
    }

    /// Returns a pointer to VesSDiffdef specified by local index.
    ///
    /// \param lidx Local index of the surface diffusion.
    VesSDiffdef& vessurfdiffdef(vessdiff_local_id lidx) const;

    /// Return the local surface diffusion index for global index argument.
    ///
    /// \param gidx Global index of the surface diffusion.
    inline vessdiff_local_id vessurfdiffG2L(vessdiff_global_id gidx) const noexcept {
        return pVesSDiff_G2L[gidx];
    }

    /// Return the global surface diffusion index for local index argument.
    ///
    /// \param lidx Local index of the surface diffusion.
    inline vessdiff_global_id vessurfdiffL2G(vessdiff_local_id lidx) const noexcept {
        return pVesSDiff_L2G[lidx];
    }

    ////////////////////////////////////////////////////////////////////////

    ///// For convenience. Why? Less convoluted way to get the required
    /// information??
    const Statedef& statedef() const noexcept {
        return pStatedef;
    }

  private:
    const Statedef& pStatedef;
    const vesicle_global_id pIdx;
    const std::string pName;
    const double pDiameter;

    // The diffusion constant
    const double pDcst;

    // Keep track of whether setup methods have been called.
    bool pSetupRefsdone{false};
    bool pSetupIndsdone{false};

    // The enclosed surface systems, stored as strings
    util::flat_set<std::string> pVssys;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Number of species embedded in vesicle surface (_V)
    ///
    uint pSpecsN_V{0};

    // Note, pSpecsN_S, _I and _O _L are not needed because we use global indices
    // for those

    // Multiple 'outer' volumes can be supported so just store the
    // global number of species
    const uint pSpecsN_global;

    /// Table to resolve species index (global -> local).
    util::strongid_vector<spec_global_id, spec_local_id> pSpec_G2L;

    /// Table to resolve species index (local -> global).
    util::strongid_vector<spec_local_id, spec_global_id> pSpec_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: LINK SPECIES
    ////////////////////////////////////////////////////////////////////////

    // Multiple 'outer' volumes can be supported so just store the
    // global number of species
    const uint pLinkSpecsN_global;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VESICLE SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of surface reactions occurring in patch.
    uint pVesSReacsN{0};

    /// Table to resolve reaction rule indices (global -> local).
    util::strongid_vector<vessreac_global_id, vessreac_local_id> pVesSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    util::strongid_vector<vessreac_local_id, vessreac_global_id> pVesSReac_L2G;

    inline uint _IDX_VesSReac_V_Spec(vessreac_local_id srlidx) const noexcept {
        return countSpecs() * srlidx.get();
    }
    inline uint _IDX_VesSReac_V_Spec(vessreac_local_id srlidx,
                                     spec_local_id splidx) const noexcept {
        return (countSpecs() * srlidx.get()) + splidx.get();
    }
    inline uint _IDX_VesSReac_L_Spec_global(vessreac_local_id srlidx) const noexcept {
        return countLinkSpecs_V() * srlidx.get();
    }
    inline uint _IDX_VesSReac_L_Spec_global(vessreac_local_id srlidx,
                                            linkspec_global_id spgidx) const noexcept {
        return (countLinkSpecs_V() * srlidx.get()) + spgidx.get();
    }
    inline uint _IDX_VesSReac_O_Spec_global(vessreac_local_id srlidx) const noexcept {
        return countSpecs_O() * srlidx.get();
    }
    inline uint _IDX_VesSReac_O_Spec_global(vessreac_local_id srlidx,
                                            spec_global_id spgidx) const noexcept {
        return (countSpecs_O() * srlidx.get()) + spgidx.get();
    }
    inline uint _IDX_VesSReac_S_Spec_global(vessreac_local_id srlidx) const noexcept {
        return countSpecs_S() * srlidx.get();
    }
    inline uint _IDX_VesSReac_S_Spec_global(vessreac_local_id srlidx,
                                            spec_global_id spgidx) const noexcept {
        return (countSpecs_S() * srlidx.get()) + spgidx.get();
    }
    inline uint _IDX_VesSReac_I_Spec_global(vessreac_local_id srlidx) const noexcept {
        return countSpecs_I() * srlidx.get();
    }
    inline uint _IDX_VesSReac_I_Spec_global(vessreac_local_id srlidx,
                                            spec_global_id spgidx) const noexcept {
        return (countSpecs_I() * srlidx.get()) + spgidx.get();
    }

    inline uint _IDX_VesSReac_Vdep_Spec_global(vessreac_local_id srlidx,
                                               spec_global_id spgidx) const noexcept {
        return (countSpecs_VDep() * srlidx.get()) + spgidx.get();
    }

    std::vector<int> pVesSReac_DEP_V_Spec;
    std::vector<int> pVesSReac_DEP_L_Spec;
    std::vector<int> pVesSReac_DEP_O_Spec;
    std::vector<int> pVesSReac_DEP_S_Spec;
    std::vector<uint> pVesSReac_LHS_V_Spec;
    std::vector<uint> pVesSReac_LHS_L_Spec;
    std::vector<uint> pVesSReac_LHS_O_Spec;
    std::vector<uint> pVesSReac_LHS_S_Spec;
    std::vector<int> pVesSReac_UPD_V_Spec;
    std::vector<int> pVesSReac_UPD_L_Spec;
    std::vector<int> pVesSReac_UPD_O_Spec;
    std::vector<int> pVesSReac_UPD_S_Spec;
    std::vector<uint> pVesSReac_UPD_I_Spec;

    std::vector<uint> pVesSReac_VDEP_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: EXOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of exocytotic reactions occurring in patch.
    uint pExocytosisN{};

    /// Table to resolve exocytotic reaction rule indices (global -> local).
    util::strongid_vector<exocytosis_global_id, exocytosis_local_id> pExocytosis_G2L;

    /// Table to resolve exocytotic reaction rule indices (local -> global).
    util::strongid_vector<exocytosis_local_id, exocytosis_global_id> pExocytosis_L2G;

    // Table of the K-constants of the exocytotic reac rules in this patch
    util::strongid_vector<exocytosis_local_id, double> pExocytosisKcst;

    // Table of 'active' flags on the exocytotic reaction rules.
    util::strongid_vector<exocytosis_local_id, uint> pExocytosisFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of surface diffusion rules occurring in vesicle.
    uint pVesSDiffsN{};

    // Table to resolve diffusion rule indices (global -> local).
    util::strongid_vector<vessdiff_global_id, vessdiff_local_id> pVesSDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
    util::strongid_vector<vessdiff_local_id, vessdiff_global_id> pVesSDiff_L2G;

    inline uint _IDX_VesSDiff_Spec(vessdiff_local_id sdiff, spec_local_id splidx) const noexcept {
        return (countSpecs() * sdiff.get()) + splidx.get();
    }

    std::vector<uint> pVesSDiff_DEP_Spec;
    util::strongid_vector<vessdiff_local_id, spec_local_id> pVesSDiff_LIG;
};

}  // namespace steps::solver
