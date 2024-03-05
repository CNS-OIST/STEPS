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

#include <gsl>

#include "fwd.hpp"
#include "geom/fwd.hpp"
#include "solver/complexeventsdef.hpp"
#include "util/collections.hpp"

namespace steps::solver {

/// Compdef object defines a compartment object.
class Compdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the compartment.
    /// \param c Associated Comp object.
    Compdef(Statedef& sd, comp_global_id idx, wm::Comp& c);

    Compdef(const Compdef&) = delete;
    Compdef& operator=(const Compdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Return the volume of this compartment.
    inline double vol() const noexcept {
        return pVol;
    }

    /// Return the global index of this compartment
    inline comp_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of this compartment.
    inline const std::string& name() const noexcept {
        return pName;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup all reference.
    void setup_references();

    /// Setup all indices.
    void setup_indices();

    /// Add a species with global index gidx.
    ///
    /// \param gidx Index of the species.
    void addSpec(spec_global_id gidx);

    /// Add an inner patch.
    ///
    /// \param p Pointer to the inner patch.
    void addIPatchdef(Patchdef& p);

    /// Add an outer patch.
    ///
    /// \param p Pointer to the outer patch.
    void addOPatchdef(Patchdef& p);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: COMPARTMENT
    ////////////////////////////////////////////////////////////////////////

    /// Set the volume of the compartment.
    ///
    /// \param v Volume of the compartment.
    void setVol(double v);

    /// Reset count, flags members of this compartment. Called when reset()
    /// method in solver object is executed.
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of species that can occur in this compartment
    inline uint countSpecs() const noexcept {
        return pSpecsN;
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

    /// Returns pointer to flags on species for this compartment.
    inline const auto& flags() const noexcept {
        return pPoolFlags;
    }

    static const uint CLAMPED = 1;

    /// Return whether a species, specified by local index argument, is
    /// clamped or not.
    ///
    /// \param sildx Local index of the species.
    inline bool clamped(spec_local_id slidx) const noexcept {
        return pPoolFlags[slidx] & CLAMPED;
    }

    /// Return pointer to species' counts in this compartment.
    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEXES
    ////////////////////////////////////////////////////////////////////////

    inline uint countComplexes() const noexcept {
        return pComplexStates.size();
    }

    inline const std::unordered_map<complex_individual_id, ComplexState>& complexStates(
        complex_global_id idx) const {
        return pComplexStates[idx];
    }

    void updateComplex(complex_global_id cmplIdx,
                       complex_individual_id stIdx,
                       const util::strongid_vector<complex_substate_id, int>& upd);

    void removeComplex(complex_global_id cmplIdx, complex_individual_id stIdx);

    void addComplex(complex_global_id cmplIdx,
                    complex_individual_id stIdx,
                    const util::strongid_vector<complex_substate_id, uint>& init);

    void clearComplexFilterUpdates();

    std::shared_ptr<ComplexFilter> GetFilter(
        const complex_global_id& cmplIdx,
        const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>&
            filts);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Set the species count of species specified by local index argument.
    ///
    /// \param slidx Local index of the species.
    /// \param count Count of molecules of a species.
    void setCount(spec_local_id slidx, double count);

    /// Clamp or unclamp species specified by local index argument
    ///
    /// \param slidx Local index of the species.
    /// \param clamp Flag to set if clamping or not.
    void setClamped(spec_local_id slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of reactions that can occur in this compartment.
    inline uint countReacs() const noexcept {
        return pReacsN;
    }

    /// Return the local reaction index for global index argument.
    ///
    /// \param gidx Global index of the reaction.
    reac_local_id reacG2L(reac_global_id gidx) const noexcept {
        return pReac_G2L[gidx];
    }

    /// Return the global reaction index for local index argument.
    ///
    /// \param lidx Local index of the reaction.
    reac_global_id reacL2G(reac_local_id lidx) const noexcept {
        return pReac_L2G[lidx];
    }

  private:
    template <typename T, typename ReacId>
    static constexpr gsl::span<const T> to_span(const std::vector<T>& container,
                                                ReacId lidx,
                                                uint size) {
        const auto begin = container.data() + lidx.get() * size;
        return {begin, begin + size};
    }

    template <typename StrongId, typename T, typename ReacId>
    static constexpr auto to_strong_span(const std::vector<T>& container, ReacId lidx, uint size) {
        return util::make_strong_random_accessor<StrongId>(to_span(container, lidx, size));
    }

  public:
    /**
     * \param rlidx Local index of the reaction
     * \return lhs array of reaction specified by local index argument
     */
    strongid_span<spec_local_id, const uint> reac_lhs(reac_local_id rlidx) const;

    /**
     * \param rlidx Local index of the reaction.
     * \return update array of reaction specified by local index argument.
     */
    strongid_span<spec_local_id, const int> reac_upd(reac_local_id rlidx) const;

    /// Return the local index of species of reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    int reac_dep(reac_local_id rlidx, spec_local_id slidx) const;

    /// Return a reference to reaction definition object (type Reacdef)
    /// specified by local index.
    ///
    /// \param rlidx Local index of the reaction.
    Reacdef& reacdef(reac_local_id rlidx) const;

    /// Return pointer to flags on reactions for this compartment.
    inline const auto& rflags() const noexcept {
        return pReacFlags;
    }

    static const uint INACTIVATED = 1;

    /// Return whether a reaction, specified by local index argument, is
    /// active or not.
    ///
    /// \param rlidx Local index of the reaction.
    inline bool active(reac_local_id rlidx) const noexcept {
        return !(pReacFlags[rlidx] & INACTIVATED);
    }

    /// Return the kcst for a reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    inline double kcst(reac_local_id rlidx) const noexcept {
        return pReacKcst[rlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of complex reactions that can occur in this compartment.
    inline uint countComplexReacs() const noexcept {
        return pComplexReacsN;
    }

    /// Return a pointer to complex reaction definition object (type Reacdef)
    /// specified by local index.
    ///
    /// \param rlidx Local index of the reaction.
    ComplexReacdef& complexreacdef(complexreac_local_id rlidx) const;

    /// Return the local complex reaction index for global index argument.
    ///
    /// \param gidx Global index of the complex reaction.
    complexreac_local_id complexreacG2L(complexreac_global_id gidx) const noexcept {
        return pComplexReac_G2L[gidx];
    }

    /// \param rlidx Local index of the complex reaction
    /// \return lhs array of complex reaction specified by local index argument
    strongid_span<spec_local_id, const uint> complexreac_lhs(complexreac_local_id rlidx) const;

    /// \param rlidx Local index of the complex reaction.
    /// \return update array of complex reaction specified by local index argument.
    strongid_span<spec_local_id, const int> complexreac_upd(complexreac_local_id rlidx) const;

    /// Return the local index of species of complex reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    int complexreac_dep(complexreac_local_id rlidx, spec_local_id slidx) const;

    /// Return whether a reaction, specified by local index argument, is
    /// active or not.
    ///
    /// \param rlidx Local index of the reaction.
    inline bool active(complexreac_local_id rlidx) const noexcept {
        return !(pComplexReacFlags[rlidx] & INACTIVATED);
    }

    /// Return the kcst for a reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    inline double kcst(complexreac_local_id rlidx) const noexcept {
        return pComplexReacKcst[rlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the reaction kcst for reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    /// \param kcst Rate constant of the reaction.
    void setKcst(reac_local_id rlidx, double kcst);

    /// Activate or inactivate a reaction specified by local index argument.
    ///
    /// \param rlidx Local index of the reaction.
    /// \param Flag to activate or inactivate a reaction.
    void setActive(reac_local_id rlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the complex reaction kcst for reaction specified by local index
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param kcst Rate constant of the reaction.
    void setComplexReacKcst(complexreac_local_id rlidx, double kcst);

    /// Activate or inactivate a complex reaction specified by local index argument.
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param Flag to activate or inactivate a complex reaction.
    void setComplexReacActive(complexreac_local_id rlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of diffusion rules for this compartment.
    inline uint countDiffs() const noexcept {
        return pDiffsN;
    }

    /// Returns a reference to Diffdef specified by local index.
    ///
    /// \param dlidx Local index of the difusion.
    Diffdef& diffdef(diff_local_id dlidx) const;

    /// Return the local diffusion index for global index argument.
    ///
    /// \param gidx Global index of the diffusion.
    diff_local_id diffG2L(diff_global_id gidx) const noexcept {
        return pDiff_G2L[gidx];
    }

    /// Return the global diffusion index for local index argument.
    ///
    /// \param lidx Local index of the diffusion.
    diff_global_id diffL2G(diff_local_id lidx) const noexcept {
        return pDiff_L2G[lidx];
    }

    /// Return the local index of species of diffusion specified by
    /// local index argument.
    ///
    /// \param dlidx Local index of the diffusion rule.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    uint diff_dep(diff_local_id dlidx, spec_local_id slidx) const;

    /// Return the rate constant of diffusion by local index argument.
    ///
    /// \param dlidx Local index of the diffusion.
    inline double dcst(diff_local_id dlidx) const noexcept {
        return pDiffDcst[dlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Set the diffusion dcst for diffusion rule specified by local index
    ///
    /// \param dlidx Local index of the diffusion.
    /// \param dcst Rate constant of the diffusion.
    void setDcst(diff_local_id dlidx, double dcst);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE BINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicle binding rules for this compartment.
    inline uint countVesBinds() const noexcept {
        return pVesBindsN;
    }

    /// Return the local ves bind index for global index argument.
    ///
    /// \param gidx Global index of the ves bind.
    vesbind_local_id vesbindG2L(vesbind_global_id gidx) const noexcept {
        return pVesBind_G2L[gidx];
    }

    /// Return the global ves bind index for local index argument.
    ///
    /// \param lidx Local index of the ves bind.
    vesbind_global_id vesbindL2G(vesbind_local_id lidx) const noexcept {
        return pVesBind_L2G[lidx];
    }

    /// Returns a reference to VesBinddef specified by local index.
    ///
    /// \param vblidx Local index of the vesicle binding.
    VesBinddef& vesbinddef(vesbind_local_id vblidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VESICLE UNBINDINGS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of vesicle unbinding rules for this compartment.
    inline uint countVesUnbinds() const noexcept {
        return pVesUnbindsN;
    }

    /// Return the local ves unbind index for global index argument.
    ///
    /// \param gidx Global index of the ves unbind.
    vesunbind_local_id vesunbindG2L(vesunbind_global_id gidx) const noexcept {
        return pVesUnbind_G2L[gidx];
    }

    /// Return the global vesunbind index for local index argument.
    ///
    /// \param lidx Local index of the vesunbind.
    vesunbind_global_id vesunbindL2G(vesunbind_local_id lidx) const noexcept {
        return pVesUnbind_L2G[lidx];
    }

    /// Returns a reference to VesUnbinddef specified by local index.
    ///
    /// \param dlidx Local index of the vesicle unbinding.
    VesUnbinddef& vesunbinddef(vesunbind_local_id local_id) const;

    ////////////////////////////////////////////////////////////////////////

    inline const Statedef& statedef() const noexcept {
        return pStatedef;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

    // A pointer to the state definition.
    Statedef& pStatedef;

    // The string identifier of the compartment
    const std::string pName;

    // The volume of the compartment
    double pVol;

    // The global index of the compartment.
    const comp_global_id pIdx;

    // The enclosed volume systems, stored as strings
    const util::flat_set<std::string>& pCvsys;

    // Keep track of whether setup_ has been called.
    bool pSetupRefsdone{false};
    bool pSetupIndsdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    // Number of species that can appear in this compartment.
    uint pSpecsN{};

    // Table to resolve species index (global -> local).
    util::strongid_vector<spec_global_id, spec_local_id> pSpec_G2L;

    // Table to resolve species index (local -> global).
    util::strongid_vector<spec_local_id, spec_global_id> pSpec_L2G;

    // Table of the populations of the species in this compartment.
    util::strongid_vector<spec_local_id, double> pPoolCount;

    // Table of 'clamped' flags on the species
    util::strongid_vector<spec_local_id, uint> pPoolFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: COMPLEXES
    ////////////////////////////////////////////////////////////////////////

    // Table of the complexes populations.
    util::strongid_vector<complex_global_id,
                          std::unordered_map<complex_individual_id, ComplexState>>
        pComplexStates;

    util::strongid_vector<complex_global_id,
                          util::strongid_vector<complex_filter_id, std::shared_ptr<ComplexFilter>>>
        pFilters;

    util::strongid_vector<
        complex_global_id,
        std::unordered_map<
            std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>,
            std::shared_ptr<ComplexFilter>,
            FilterHash>>
        pFiltersMap;

    ////////////////////////////////////////////////////////////////////////
    // DATA: REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of reaction rules occurring in this compartment.
    uint pReacsN{};

    // Table to resolve reaction rule indices (global -> local).
    util::strongid_vector<reac_global_id, reac_local_id> pReac_G2L;

    // Table to resolve reaction rule indices (local -> global).
    util::strongid_vector<reac_local_id, reac_global_id> pReac_L2G;

    // Table of the K-constants of the reaction rules in this compartment
    util::strongid_vector<reac_local_id, double> pReacKcst;

    // Table of 'active' flags on the reaction rules.
    util::strongid_vector<reac_local_id, uint> pReacFlags;

    inline uint _IDX_Reac_Spec(reac_local_id reac, spec_local_id spec) const {
        return (countSpecs() * reac.get()) + spec.get();
    }

    std::vector<int> pReac_DEP_Spec;
    std::vector<uint> pReac_LHS_Spec;
    std::vector<int> pReac_UPD_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: COMPLEX REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Number of complex reaction rules occurring in this compartment.
    uint pComplexReacsN{};

    /// Table to resolve complex reaction rule indices (global -> local).
    util::strongid_vector<complexreac_global_id, complexreac_local_id> pComplexReac_G2L;

    /// Table to resolve complex reaction rule indices (local -> global).
    util::strongid_vector<complexreac_local_id, complexreac_global_id> pComplexReac_L2G;

    /// Table of the K-constants of the complex reaction rules in this compartment
    util::strongid_vector<complexreac_local_id, double> pComplexReacKcst;

    /// Table of 'active' flags on the complex reaction rules.
    util::strongid_vector<complexreac_local_id, uint> pComplexReacFlags;

    inline uint _IDX_ComplexReac_Spec(complexreac_local_id reac, spec_local_id spec) const {
        return (countSpecs() * reac.get()) + spec.get();
    }

    std::vector<int> pComplexReac_DEP_Spec;
    std::vector<uint> pComplexReac_LHS_Spec;
    std::vector<int> pComplexReac_UPD_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of diffusion rules occuring in compartment.
    uint pDiffsN{};

    // Table to resolve diffusion rule indices (global -> local).
    util::strongid_vector<diff_global_id, diff_local_id> pDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
    util::strongid_vector<diff_local_id, diff_global_id> pDiff_L2G;

    // Table of the D-constants of the diffusion rules in this compartment
    util::strongid_vector<diff_local_id, double> pDiffDcst;

    inline uint _IDX_Diff_Spec(diff_local_id diff, spec_local_id spec) const noexcept {
        return (countSpecs() * diff.get()) + spec.get();
    }

    std::vector<uint> pDiff_DEP_Spec;
    util::strongid_vector<diff_local_id, spec_local_id> pDiff_LIG;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VESICLE BINDING RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of vesicle binding rules occuring in compartment.
    uint pVesBindsN{};

    // Table to resolve vesicle binding rule indices (global -> local).
    util::strongid_vector<vesbind_global_id, vesbind_local_id> pVesBind_G2L;

    // Table to resolve vesicle binding rule indices (local -> global).
    util::strongid_vector<vesbind_local_id, vesbind_global_id> pVesBind_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VESICLE UNBINDING RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of vesicle unbinding rules occuring in compartment.
    uint pVesUnbindsN{};

    // Table to resolve vesicle unbinding rule indices (global -> local).
    util::strongid_vector<vesunbind_global_id, vesunbind_local_id> pVesUnbind_G2L;

    // Table to resolve vesicle unbinding rule indices (local -> global).
    util::strongid_vector<vesunbind_local_id, vesunbind_global_id> pVesUnbind_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    // A list of patches surrounding
    std::vector<Patchdef*> pIPatches;
    std::vector<Patchdef*> pOPatches;
};

}  // namespace steps::solver
