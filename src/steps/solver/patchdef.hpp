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

#include "complexeventsdef.hpp"
#include "endocyticzonedef.hpp"
#include "fwd.hpp"
#include "geom/fwd.hpp"
#include "model/complexevents.hpp"

namespace steps::solver {

/// Defined patch object.
class Patchdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the patch.
    /// \param p Reference to the Patch object.
    Patchdef(Statedef& sd, patch_global_id idx, wm::Patch& p);

    Patchdef(const Patchdef&) = delete;
    Patchdef& operator=(const Patchdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: PATCH
    ////////////////////////////////////////////////////////////////////////

    /// Return the area of this patch.
    inline double area() const noexcept {
        return pArea;
    }

    /// Return the global index of this patch.
    inline patch_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the patch.
    const std::string& name() const noexcept {
        return pName;
    }

    /// Return a pointer to the inner compartment.
    inline Compdef* icompdef() const noexcept {
        return pInner;
    }

    /// Return a pointer to the outer compartment.
    inline Compdef* ocompdef() const noexcept {
        return pOuter;
    }

    /// Return all endocytic zones.
    inline const auto& endocyticZones() const noexcept {
        return pEndocyticZonesdefs;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup all references.
    void setup_references();

    /// Setup all indices.
    void setup_indices();

    /// Add a species with global index gidx.
    ///
    /// \param gidx Index of the species.
    void addSpec(spec_global_id gidx);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: PATCH
    ////////////////////////////////////////////////////////////////////////

    /// Set the area of this patch.
    ///
    /// \param a Area of the patch.
    void setArea(double a);

    /// Get the area of the patch
    ///
    /// \return Area of the patch
    inline double getArea() const noexcept {
        return pArea;
    }

    /// Reset count, flags members of this patch. Called when reset()
    /// method in solver object is executed.
    ///
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of species defined for this surface patch.
    inline uint countSpecs() const noexcept {
        return pSpecsN_S;
    }

    /// Return the number of species defined for the inner compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_I() const noexcept {
        return pSpecsN_I;
    }

    /// Return the number of species defined for the outer compartment.
    /// Should not be called before Compdef::setup()
    inline uint countSpecs_O() const noexcept {
        return pSpecsN_O;
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

    /// Auxiliary function: resolves a species gidx for the inner
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    spec_local_id specG2L_I(spec_global_id gidx) const;
    /// Auxiliary function: resolves a species gidx for the outer
    /// compartment.
    ///
    /// \param gidx Global index of the species.
    /// \return The local index or steps::sim::shared::LIDX_UNDEFINED.
    ///
    spec_local_id specG2L_O(spec_global_id gidx) const;

    /// Return pointer to species' counts on this patch.
    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    /// Returns pointer to flags on species for this patch.
    inline const auto& flags() const noexcept {
        return pPoolFlags;
    }

    static const uint CLAMPED = 1;

    /// Return whether a species, specified by local index argument, is
    /// clamped or not.
    ///
    /// \param slidx Local index of the species.
    inline bool clamped(spec_local_id slidx) const noexcept {
        return pPoolFlags[slidx] & CLAMPED;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Set the species count of species specified by local index argument.
    ///
    /// \param slidx Local index of the species.
    /// \param count Count of species.
    void setCount(spec_local_id slidx, double count);

    /// Clamp or unclamp species specified by local index argument
    ///
    /// \param slidx Local index of the species.
    /// \param clamp Flag to clamp or unclamp species.
    void setClamped(spec_local_id slidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEXES
    ////////////////////////////////////////////////////////////////////////

    inline uint countComplexes() const noexcept {
        return pComplexStates.container().size();
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
    // DATA ACCESS: SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface reactions that can occur on this
    /// patch.
    inline uint countSReacs() const noexcept {
        return pSReacsN;
    }

    /// Return the local surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the surface reaction.
    inline sreac_local_id sreacG2L(sreac_global_id gidx) const noexcept {
        return pSReac_G2L[gidx];
    }

    /// Return the global surface reaction index for local index argument.
    ///
    /// \param lidx Local index of the surface reaction.
    inline sreac_global_id sreacL2G(sreac_local_id lidx) const noexcept {
        return pSReac_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type SReacdef)
    /// specified by local index.
    SReacdef& sreacdef(sreac_local_id lidx) const;

    /// Warning: these methods perform no error checking!
    int sreac_dep_I(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pSReac_DEP_I_Spec[splidx.get() + (srlidx.get() * countSpecs_I())];
    }

    int sreac_dep_S(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pSReac_DEP_S_Spec[splidx.get() + (srlidx.get() * countSpecs())];
    }

    int sreac_dep_O(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pSReac_DEP_O_Spec[splidx.get() + (srlidx.get() * countSpecs_O())];
    }

  private:
    /// internal utility function
    template <typename T, typename ReacLocalId>
    static constexpr gsl::span<const T> to_span(const std::vector<T>& container,
                                                ReacLocalId lidx,
                                                uint size) {
        const auto begin = container.data() + lidx.get() * size;
        return {begin, begin + size};
    }

    template <typename StrongId, typename T, typename ReacLocalId>
    static constexpr auto to_strong_span(const std::vector<T>& container,
                                         ReacLocalId lidx,
                                         uint size) {
        return util::make_strong_random_accessor<StrongId>(to_span(container, lidx, size));
    }

  public:
    ///
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const uint> sreac_lhs_I(sreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pSReac_LHS_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const uint> sreac_lhs_S(sreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pSReac_LHS_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const uint> sreac_lhs_O(sreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pSReac_LHS_O_Spec, lidx, countSpecs_O());
    }

    ///
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const int> sreac_upd_I(sreac_local_id lidx) const noexcept {
        return to_span(pSReac_UPD_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const int> sreac_upd_S(sreac_local_id lidx) const noexcept {
        return to_span(pSReac_UPD_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const int> sreac_upd_O(sreac_local_id lidx) const noexcept {
        return to_span(pSReac_UPD_O_Spec, lidx, countSpecs_O());
    }

    /// Return pointer to flags on surface reactions for this patch.
    inline const auto& srflags() const noexcept {
        return pSReacFlags;
    }

    static const uint INACTIVATED = 1;

    /// Return whether a surface reaction, specified by local index argument,
    /// is active or not.
    ///
    /// \param rlidx Local index of the surface reaction.
    inline bool active(sreac_local_id rlidx) const noexcept {
        return !(pSReacFlags[rlidx] & INACTIVATED);
    }

    /// Return the kcst for a surface reaction specified by local index
    ///
    /// \param rlidx Local index of the surface reaction.
    inline double kcst(sreac_local_id rlidx) const noexcept {
        return pSReacKcst[rlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of complex surface reactions that can occur in this compartment.
    inline uint countComplexSReacs() const noexcept {
        return pComplexSReacsN;
    }

    /// Return a pointer to complex surface reaction definition object
    /// specified by local index.
    ///
    /// \param rlidx Local index of the reaction.
    ComplexSReacdef& complexsreacdef(complexsreac_local_id rlidx) const;

    /// Return the local complex surface reaction index for global index argument.
    ///
    /// \param gidx Global index of the complex reaction.
    complexsreac_local_id complexsreacG2L(complexsreac_global_id gidx) const noexcept {
        return pComplexSReac_G2L[gidx];
    }

    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const uint> complexsreac_lhs_I(
        complexsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pComplexSReac_LHS_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const uint> complexsreac_lhs_S(
        complexsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pComplexSReac_LHS_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const uint> complexsreac_lhs_O(
        complexsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pComplexSReac_LHS_O_Spec, lidx, countSpecs_O());
    }

    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const int> complexsreac_upd_I(
        complexsreac_local_id lidx) const noexcept {
        return to_span(pComplexSReac_UPD_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const int> complexsreac_upd_S(
        complexsreac_local_id lidx) const noexcept {
        return to_span(pComplexSReac_UPD_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const int> complexsreac_upd_O(
        complexsreac_local_id lidx) const noexcept {
        return to_span(pComplexSReac_UPD_O_Spec, lidx, countSpecs_O());
    }

    /// Return the local index of species of complex reaction specified by
    /// local index argument.
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param slidx Local index of the species.
    int complexsreac_dep_I(complexsreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pComplexSReac_DEP_I_Spec[splidx.get() + (srlidx.get() * countSpecs_I())];
    }

    int complexsreac_dep_S(complexsreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pComplexSReac_DEP_S_Spec[splidx.get() + (srlidx.get() * countSpecs())];
    }

    int complexsreac_dep_O(complexsreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return pComplexSReac_DEP_O_Spec[splidx.get() + (srlidx.get() * countSpecs_O())];
    }

    /// Return whether a reaction, specified by local index argument, is
    /// active or not.
    ///
    /// \param rlidx Local index of the reaction.
    inline bool active(complexsreac_local_id rlidx) const noexcept {
        return !(pComplexSReacFlags[rlidx] & INACTIVATED);
    }

    /// Return the kcst for a reaction specified by local index
    ///
    /// \param rlidx Local index of the reaction.
    inline double kcst(complexsreac_local_id rlidx) const noexcept {
        return pComplexSReacKcst[rlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: ENDOCYTOTIC REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of endocytotic reactions that can occur on this
    /// patch.
    inline uint countEndocytosis() const noexcept {
        return pEndocytosisN;
    }

    /// Return the global endocytotic reaction index for global index argument.
    ///
    /// \param gidx Global index of the endocytotic reaction.
    inline endocytosis_local_id endocytosisG2L(endocytosis_global_id gidx) const noexcept {
        return pEndocytosis_G2L[gidx];
    }

    /// Return the local endocytotic reaction index for global index argument.
    ///
    /// \param lidx Local index of the endocytotic reaction.
    inline endocytosis_global_id endocytosisL2G(endocytosis_local_id lidx) const noexcept {
        return pEndocytosis_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type Endocytosisdef)
    /// specified by local index.
    Endocytosisdef& endocytosisdef(endocytosis_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: RAFT GENESIS RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of raft geneses reactions that can occur on this
    /// patch.
    inline uint countRaftGens() const noexcept {
        return pRaftGenN;
    }

    /// Return the local raft genesis index for global index argument.
    ///
    /// \param gidx Global index of the raft genesis
    inline raftgen_local_id raftgenG2L(raftgen_global_id gidx) const noexcept {
        return pRaftGen_G2L[gidx];
    }

    /// Return the global raft genesis index for local index argument.
    ///
    /// \param lidx Local index of the raft genesis
    inline raftgen_global_id raftgenL2G(raftgen_local_id lidx) const noexcept {
        return pRaftGen_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type RaftGendef)
    /// specified by local index.
    RaftGendef& raftgendef(raftgen_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SURFACE DIFFUSION
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of surface diffusion rules for this patch.
    inline uint countSurfDiffs() const noexcept {
        return pSurfDiffsN;
    }

    /// Returns a reference to Diffdef specified by local index.
    ///
    /// \param dlidx Local index of the surface diffusion.
    SurfDiffdef& surfdiffdef(surfdiff_local_id dlidx) const;

    /// Return the local surface diffusion index for global index argument.
    ///
    /// \param gidx Global index of the surface diffusion.
    surfdiff_local_id surfdiffG2L(surfdiff_global_id gidx) const noexcept {
        return pSurfDiff_G2L[gidx];
    }

    /// Return the global surface diffusion index for local index argument.
    ///
    /// \param lidx Local index of the surface diffusion.
    surfdiff_global_id surfdiffL2G(surfdiff_local_id lidx) const noexcept {
        return pSurfDiff_L2G[lidx];
    }

    /// Return the local index of species of surface diffusion specified by
    /// local index argument.
    ///
    /// \param dlidx Local index of the surface diffusion rule.
    /// \param slidx Local index of the species.
    /// \todo make sure this is correct.
    uint surfdiff_dep(surfdiff_local_id dlidx, spec_local_id slidx) const;

    /// Return the rate constant of surface diffusion by local index argument.
    ///
    /// \param dlidx Local index of the surface diffusion.
    inline double dcst(surfdiff_local_id dlidx) const noexcept {
        return pSurfDiffDcst[dlidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of ohmic currents defined in this patch.
    inline uint countOhmicCurrs() const noexcept {
        return pOhmicCurrsN;
    }

    /// Return the local index of ohmic current with global index gidx.
    inline ohmiccurr_local_id ohmiccurrG2L(ohmiccurr_global_id gidx) const noexcept {
        return pOhmicCurr_G2L[gidx];
    }

    /// Return the global index of ohmic current with local index lidx.
    inline ohmiccurr_global_id ohmiccurrL2G(ohmiccurr_local_id lidx) const noexcept {
        return pOhmicCurr_L2G[lidx];
    }

    /// Return a reference to ohmic current definition object (type OhmicCurrdef)
    /// specified by local index oclidx.
    OhmicCurrdef& ohmiccurrdef(ohmiccurr_local_id oclidx) const;

    /// Return dependency information of ohmic current with local index oclidx
    /// on species with local index splidx.
    int ohmiccurr_dep_S(ohmiccurr_local_id oclidx, spec_local_id splidx) const;

    /// Return the local index of the channel state associated with the
    /// ohmic current with local index oclidx
    spec_local_id ohmiccurr_chanstate(ohmiccurr_local_id oclidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    /// Return the number of GHK currents defined in this patch.
    inline uint countGHKcurrs() const noexcept {
        return pGHKcurrsN;
    }

    /// Return the local index of GHK current with global index gidx.
    inline ghkcurr_local_id ghkcurrG2L(ghkcurr_global_id gidx) const noexcept {
        return pGHKcurr_G2L[gidx];
    }

    /// Return the global index of GHK current with local index lidx.
    inline ghkcurr_global_id ghkcurrL2G(ghkcurr_local_id lidx) const noexcept {
        return pGHKcurr_L2G[lidx];
    }

    /// Return a reference to GHK current definition object (type GHKcurrdef)
    /// specified by local index lidx.
    GHKcurrdef& ghkcurrdef(ghkcurr_local_id ghklidx) const;

    /// Return dependency information of ghk current with local index ghklidx
    /// on species with local index splidx.
    int ghkcurr_dep_S(ghkcurr_local_id ghklidx, spec_local_id splidx) const;

    /// Return the local index of the channel state associated with the
    /// ghk current with local index ghklidx
    spec_local_id ghkcurr_chanstate(ghkcurr_local_id ghklidx) const;

    /// Return the local index of the ION associated with the ghk current
    // with local index ghklidx
    spec_local_id ghkcurr_ion(ghkcurr_local_id ghklidx) const;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: VOLTAGE-DEPENDENT REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Return the total number of voltage-dependent reactions that can occur on
    /// this patch.
    inline uint countVDepSReacs() const noexcept {
        return pVDepSReacsN;
    }

    /// Return the local voltage-dependent reaction index for global index
    /// argument.
    ///
    /// \param gidx Global index of the voltage-dependent reaction.
    inline vdepsreac_local_id vdepsreacG2L(vdepsreac_global_id gidx) const noexcept {
        return pVDepSReac_G2L[gidx];
    }

    /// Return the global voltage-dependent reaction index for local index
    /// argument.
    ///
    /// \param lidx Local index of the voltage-dependent reaction.
    inline vdepsreac_global_id vdepsreacL2G(vdepsreac_local_id lidx) const noexcept {
        return pVDepSReac_L2G[lidx];
    }

    /// Return a reference to reaction definition object (type VDepSReacdef)
    /// specified by local index.
    VDepSReacdef& vdepsreacdef(vdepsreac_local_id lidx) const;

    /// Warning: these methods perform no error checking!
    /// \todo imcompleted.
    int vdepsreac_dep_I(vdepsreac_local_id vdsrlidx, spec_local_id splidx) const noexcept {
        return pVDepSReac_DEP_I_Spec[splidx.get() + (vdsrlidx.get() * countSpecs_I())];
    }

    int vdepsreac_dep_S(vdepsreac_local_id vdsrlidx, spec_local_id splidx) const noexcept {
        return pVDepSReac_DEP_S_Spec[splidx.get() + (vdsrlidx.get() * countSpecs())];
    }

    int vdepsreac_dep_O(vdepsreac_local_id vdsrlidx, spec_local_id splidx) const noexcept {
        return pVDepSReac_DEP_O_Spec[splidx.get() + (vdsrlidx.get() * countSpecs_O())];
    }


    /// \param lidx
    /// \return
    // Return the beginning and end of the lhs arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const uint> vdepsreac_lhs_I(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_LHS_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const uint> vdepsreac_lhs_S(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_LHS_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const uint> vdepsreac_lhs_O(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_LHS_O_Spec, lidx, countSpecs_O());
    }

    ///
    // Return the beginning and end of the update arrays of surface reaction
    // specified by local index argument.
    strongid_span<spec_local_id, const int> vdepsreac_upd_I(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_UPD_I_Spec, lidx, countSpecs_I());
    }

    strongid_span<spec_local_id, const int> vdepsreac_upd_S(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_UPD_S_Spec, lidx, countSpecs());
    }

    strongid_span<spec_local_id, const int> vdepsreac_upd_O(
        vdepsreac_local_id lidx) const noexcept {
        return to_strong_span<spec_local_id>(pVDepSReac_UPD_O_Spec, lidx, countSpecs_O());
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the kcst for surface reaction specified by local index
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param kcst Rate constant of the surface reaction.
    void setKcst(sreac_local_id srlidx, double kcst);

    /// Activate or inactivate a surface reaction specified by local index.
    ///
    /// \param srlidx Local index of the surface reaction.
    /// \param active Flag to activate / inactivate the surface reaction.
    void setActive(sreac_local_id srlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: COMPLEX REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the complex surface reaction kcst for reaction specified by local index
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param kcst Rate constant of the reaction.
    void setComplexSReacKcst(complexsreac_local_id rlidx, double kcst);

    /// Activate or inactivate a complex surface reaction specified by local index argument.
    ///
    /// \param rlidx Local index of the complex reaction.
    /// \param Flag to activate or inactivate a complex reaction.
    void setComplexSReacActive(complexsreac_local_id rlidx, bool active);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: ENDOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Set the kcst for endocytotic reaction specified by local index
    ///
    /// \param endolidx Local index of the endocytotic reaction.
    /// \param kcst Rate constant of the endocytotic reaction.
    void setEndoKcst(endocytosis_local_id endolidx, double kcst);

    /// Activate or inactivate a endocytotic reaction specified by local index.
    ///
    /// \param endolidx Local index of the endocytotic reaction.
    /// \param active Flag to activate / inactivate the endocytotic reaction.
    void setEndoActive(endocytosis_local_id endolidx, bool active);

    inline const Statedef& statedef() const noexcept {
        return pStatedef;
    }

  private:
    ////////////////////////////////////////////////////////////////////////
    // DATA: GENERAL
    ////////////////////////////////////////////////////////////////////////

    /// A pointer to the state definition.
    const Statedef& pStatedef;

    // The string identifier of the patch
    const std::string pName;

    // The area of the patch
    double pArea;

    /// The index of the patch.
    patch_global_id pIdx;

    // geom level comps, to be used during setup and NOT later.
    const wm::Comp& pIcomp;
    const wm::Comp* pOcomp;

    const util::flat_set<std::string> pPssys;

    /// Pointer to inner compartment CompDef.
    Compdef* pInner{nullptr};

    /// Pointer to outer compartment CompDef.
    Compdef* pOuter{nullptr};

    // Keep track of whether setup methods have been called.
    bool pSetupRefsdone{false};
    bool pSetupIndsdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: SPECIES
    ////////////////////////////////////////////////////////////////////////

    /// Number of species embedded in inner volume (_I), patch (_S)
    /// and outer volume (_O).
    uint pSpecsN_I{};
    uint pSpecsN_S{};
    uint pSpecsN_O{};

    /// Table to resolve species index (global -> local).
    util::strongid_vector<spec_global_id, spec_local_id> pSpec_G2L;

    /// Table to resolve species index (local -> global).
    util::strongid_vector<spec_local_id, spec_global_id> pSpec_L2G;

    // Table of the populations of the species on this patch.
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
    // DATA: SURFACE REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of surface reactions occurring in patch.
    uint pSReacsN{};

    /// Table to resolve reaction rule indices (global -> local).
    util::strongid_vector<sreac_global_id, sreac_local_id> pSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    util::strongid_vector<sreac_local_id, sreac_global_id> pSReac_L2G;

    // Table of the K-constants of the surface reac rules in this patch
    util::strongid_vector<sreac_local_id, double> pSReacKcst;

    // Table of 'active' flags on the surface reaction rules.
    util::strongid_vector<sreac_local_id, uint> pSReacFlags;

    inline uint _IDX_SReac_I_Spec(sreac_local_id srlidx) const noexcept {
        return countSpecs_I() * srlidx.get();
    }
    inline uint _IDX_SReac_I_Spec(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return (countSpecs_I() * srlidx.get()) + splidx.get();
    }
    inline uint _IDX_SReac_S_Spec(sreac_local_id srlidx) const noexcept {
        return countSpecs() * srlidx.get();
    }
    inline uint _IDX_SReac_S_Spec(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return (countSpecs() * srlidx.get()) + splidx.get();
    }
    inline uint _IDX_SReac_O_Spec(sreac_local_id srlidx) const noexcept {
        return countSpecs_O() * srlidx.get();
    }
    inline uint _IDX_SReac_O_Spec(sreac_local_id srlidx, spec_local_id splidx) const noexcept {
        return (countSpecs_O() * srlidx.get()) + splidx.get();
    }

    std::vector<int> pSReac_DEP_I_Spec;
    std::vector<int> pSReac_DEP_S_Spec;
    std::vector<int> pSReac_DEP_O_Spec;
    std::vector<uint> pSReac_LHS_I_Spec;
    std::vector<uint> pSReac_LHS_S_Spec;
    std::vector<uint> pSReac_LHS_O_Spec;
    std::vector<int> pSReac_UPD_I_Spec;
    std::vector<int> pSReac_UPD_S_Spec;
    std::vector<int> pSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: COMPLEX SURFACE REACTION RULES
    ////////////////////////////////////////////////////////////////////////

    /// Number of complex reaction rules occuring in this compartment.
    uint pComplexSReacsN{};
    /// Table to resolve complex surface reaction rule indices (global -> local).
    util::strongid_vector<complexsreac_global_id, complexsreac_local_id> pComplexSReac_G2L;
    /// Table to resolve complex surface reaction rule indices (local -> global).
    util::strongid_vector<complexsreac_local_id, complexsreac_global_id> pComplexSReac_L2G;
    /// Table of the K-constants of the complex surface reaction rules in this compartment
    util::strongid_vector<complexsreac_local_id, double> pComplexSReacKcst;
    /// Table of 'active' flags on the complex surface reaction rules.
    util::strongid_vector<complexsreac_local_id, uint> pComplexSReacFlags;

    inline uint _IDX_ComplexSReac_I_Spec(complexsreac_local_id reac) const {
        return countSpecs_I() * reac.get();
    }
    inline uint _IDX_ComplexSReac_I_Spec(complexsreac_local_id reac, spec_local_id spec) const {
        return (countSpecs_I() * reac.get()) + spec.get();
    }
    inline uint _IDX_ComplexSReac_S_Spec(complexsreac_local_id reac) const {
        return countSpecs() * reac.get();
    }
    inline uint _IDX_ComplexSReac_S_Spec(complexsreac_local_id reac, spec_local_id spec) const {
        return (countSpecs() * reac.get()) + spec.get();
    }
    inline uint _IDX_ComplexSReac_O_Spec(complexsreac_local_id reac) const {
        return countSpecs_O() * reac.get();
    }
    inline uint _IDX_ComplexSReac_O_Spec(complexsreac_local_id reac, spec_local_id spec) const {
        return (countSpecs_O() * reac.get()) + spec.get();
    }
    std::vector<int> pComplexSReac_DEP_I_Spec;
    std::vector<int> pComplexSReac_DEP_S_Spec;
    std::vector<int> pComplexSReac_DEP_O_Spec;
    std::vector<uint> pComplexSReac_LHS_I_Spec;
    std::vector<uint> pComplexSReac_LHS_S_Spec;
    std::vector<uint> pComplexSReac_LHS_O_Spec;
    std::vector<int> pComplexSReac_UPD_I_Spec;
    std::vector<int> pComplexSReac_UPD_S_Spec;
    std::vector<int> pComplexSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: SURFACE DIFFUSION RULES
    ////////////////////////////////////////////////////////////////////////

    // Number of surface diffusion rules occurring in patch.
    uint pSurfDiffsN{};

    // Table to resolve diffusion rule indices (global -> local).
    util::strongid_vector<surfdiff_global_id, surfdiff_local_id> pSurfDiff_G2L;

    // Table to resolve diffusion rule indices (local -> global).
    util::strongid_vector<surfdiff_local_id, surfdiff_global_id> pSurfDiff_L2G;

    // Table of the D-constants of the diffusion rules in this compartment
    util::strongid_vector<surfdiff_local_id, double> pSurfDiffDcst;

    inline uint _IDX_SurfDiff_Spec(surfdiff_local_id sdlidx, spec_local_id splidx) const noexcept {
        return (pSpecsN_S * sdlidx.get()) + splidx.get();
    }
    std::vector<uint> pSurfDiff_DEP_Spec;
    util::strongid_vector<surfdiff_local_id, spec_local_id> pSurfDiff_LIG;

    ////////////////////////////////////////////////////////////////////////
    // DATA: ENDOCYTOTIC ZONES
    ////////////////////////////////////////////////////////////////////////

    std::vector<std::unique_ptr<EndocyticZonedef>> pEndocyticZonesdefs;

    ////////////////////////////////////////////////////////////////////////
    // DATA: ENDOCYTOTIC REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of endocytotic reactions occurring in patch.
    uint pEndocytosisN{};

    /// Table to resolve endocytotic reaction rule indices (global -> local).
    util::strongid_vector<endocytosis_global_id, endocytosis_local_id> pEndocytosis_G2L;

    /// Table to resolve endocytotic reaction rule indices (local -> global).
    util::strongid_vector<endocytosis_local_id, endocytosis_global_id> pEndocytosis_L2G;

    // Table of the K-constants of the endocytotic reac rules in this patch
    util::strongid_vector<endocytosis_local_id, double> pEndocytosisKcst;

    // Table of 'active' flags on the endocytotic reaction rules.
    util::strongid_vector<endocytosis_local_id, uint> pEndocytosisFlags;

    ////////////////////////////////////////////////////////////////////////
    // DATA: RAFT GENESIS
    ////////////////////////////////////////////////////////////////////////

    /// Number of raft geneses occurring in patch.
    uint pRaftGenN{};

    /// Table to resolve raft genesis rule indices (global -> local).
    util::strongid_vector<raftgen_global_id, raftgen_local_id> pRaftGen_G2L;

    /// Table to resolve raft genesis rule indices (local -> global).
    util::strongid_vector<raftgen_local_id, raftgen_global_id> pRaftGen_L2G;

    ////////////////////////////////////////////////////////////////////////
    // DATA: VOLTAGE-DEPENDENT REACTIONS
    ////////////////////////////////////////////////////////////////////////

    /// Number of voltage-dependent reactions occurring in patch.
    uint pVDepSReacsN{};

    /// Table to resolve reaction rule indices (global -> local).
    util::strongid_vector<vdepsreac_global_id, vdepsreac_local_id> pVDepSReac_G2L;

    /// Table to resolve reaction rule indices (local -> global).
    util::strongid_vector<vdepsreac_local_id, vdepsreac_global_id> pVDepSReac_L2G;

    inline uint _IDX_VDepSReac_I_Spec(vdepsreac_local_id vdsrlidx) const noexcept {
        return countSpecs_I() * vdsrlidx.get();
    }
    inline uint _IDX_VDepSReac_I_Spec(vdepsreac_local_id vdsrlidx,
                                      spec_local_id splidx) const noexcept {
        return (countSpecs_I() * vdsrlidx.get()) + splidx.get();
    }
    inline uint _IDX_VDepSReac_S_Spec(vdepsreac_local_id vdsrlidx) const noexcept {
        return countSpecs() * vdsrlidx.get();
    }
    inline uint _IDX_VDepSReac_S_Spec(vdepsreac_local_id vdsrlidx,
                                      spec_local_id splidx) const noexcept {
        return (countSpecs() * vdsrlidx.get()) + splidx.get();
    }
    inline uint _IDX_VDepSReac_O_Spec(vdepsreac_local_id vdsrlidx) const noexcept {
        return countSpecs_O() * vdsrlidx.get();
    }
    inline uint _IDX_VDepSReac_O_Spec(vdepsreac_local_id vdsrlidx,
                                      spec_local_id splidx) const noexcept {
        return (countSpecs_O() * vdsrlidx.get()) + splidx.get();
    }

    std::vector<int> pVDepSReac_DEP_I_Spec;
    std::vector<int> pVDepSReac_DEP_S_Spec;
    std::vector<int> pVDepSReac_DEP_O_Spec;
    std::vector<uint> pVDepSReac_LHS_I_Spec;
    std::vector<uint> pVDepSReac_LHS_S_Spec;
    std::vector<uint> pVDepSReac_LHS_O_Spec;
    std::vector<int> pVDepSReac_UPD_I_Spec;
    std::vector<int> pVDepSReac_UPD_S_Spec;
    std::vector<int> pVDepSReac_UPD_O_Spec;

    ////////////////////////////////////////////////////////////////////////
    // DATA: OHMIC CURRENTS
    ////////////////////////////////////////////////////////////////////////

    uint pOhmicCurrsN{};
    util::strongid_vector<ohmiccurr_global_id, ohmiccurr_local_id> pOhmicCurr_G2L;
    util::strongid_vector<ohmiccurr_local_id, ohmiccurr_global_id> pOhmicCurr_L2G;

    inline uint _IDX_OhmicCurr_Spec(ohmiccurr_local_id ohmicclidx) const noexcept {
        return countOhmicCurrs() * ohmicclidx.get();
    }

    inline uint _IDX_OhmicCurr_Spec(ohmiccurr_local_id ohmicclidx,
                                    spec_local_id splidx) const noexcept {
        return (countOhmicCurrs() * ohmicclidx.get()) + splidx.get();
    }

    std::vector<int> pOhmicCurr_DEP_Spec;
    util::strongid_vector<ohmiccurr_local_id, spec_local_id> pOhmicCurr_CHANSTATE;

    ////////////////////////////////////////////////////////////////////////
    // DATA: GHK CURRENTS
    ////////////////////////////////////////////////////////////////////////

    uint pGHKcurrsN{};
    util::strongid_vector<ghkcurr_global_id, ghkcurr_local_id> pGHKcurr_G2L;
    util::strongid_vector<ghkcurr_local_id, ghkcurr_global_id> pGHKcurr_L2G;

    inline uint _IDX_GHKcurr_Spec(ghkcurr_local_id ghkclidx) {
        return countGHKcurrs() * ghkclidx.get();
    }
    inline uint _IDX_GHKcurr_Spec(ghkcurr_local_id ghkclidx, spec_local_id splidx) const noexcept {
        return (countGHKcurrs() * ghkclidx.get()) + splidx.get();
    }

    std::vector<int> pGHKcurr_DEP_Spec;
    util::strongid_vector<ghkcurr_local_id, spec_local_id> pGHKcurr_CHANSTATE;
    util::strongid_vector<ghkcurr_local_id, spec_local_id> pGHKcurr_ION;
};

}  // namespace steps::solver
