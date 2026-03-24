#pragma once

#include <cassert>
#include <memory>
#include <numeric>
#include <optional>
#include <set>
#include <type_traits>
#include <unordered_map>
#include <vector>

#include "fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "model/fwd.hpp"
#include "model/ghkcurr.hpp"
#include "model/spec.hpp"
#include "model/vdepsreac.hpp"
#include "mpi/dist/tetopsplit/definition/efield.hpp"
#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "solver/vdepsreacdef.hpp"
#include "sreacdef.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {


/**
 * \brief State definition of a patch.
 *
 * The Patchdef class defines the sub biochemical container of a patch.
 * It provides the global and local indexing for species, reactions
 * and diffusions in the compartment.
 *
 */
class Patchdef {
  public:
    using AllModSReacID = std::variant<model::surface_reaction_id,
                                       model::complex_surface_reaction_id,
                                       model::vdep_surface_reaction_id,
                                       model::vdep_complex_surface_reaction_id,
                                       model::ghk_current_id,
                                       model::complex_ghk_current_id>;
    using AllContSReacID = std::variant<container::surface_reaction_id,
                                        container::complex_surface_reaction_id,
                                        container::vdep_surface_reaction_id,
                                        container::vdep_complex_surface_reaction_id,
                                        container::ghk_surface_reaction_id,
                                        container::complex_ghk_surface_reaction_id>;
    using AllModCurrentID = std::variant<model::ohmic_current_id,
                                         model::complex_ohmic_current_id,
                                         model::ghk_current_id,
                                         model::complex_ghk_current_id>;
    using AllContCurrentID = std::variant<container::ohmic_current_id,
                                          container::complex_ohmic_current_id,
                                          container::ghk_current_id,
                                          container::complex_ghk_current_id>;

    Patchdef(const Statedef& statedef,
             const DistPatch& patch,
             container::patch_id container_patch_id);

    /// Reset the def-object values to model defaults
    void reset();

    inline const DistPatch& patch() const noexcept {
        return pPatch_;
    }

    inline const model::patch_id& getID() const noexcept {
        return model_patch_;
    }

    inline container::patch_id getIdx() const noexcept {
        return container_patch_id_;
    }

    inline model::compartment_id getInnerCompId() const noexcept {
        return inner_compartment_id_;
    }

    template <typename MReacID>
    typename ModID2ContID<MReacID>::type getReacIdx(MReacID reac) const {
        using ReacID = typename ModID2ContID<MReacID>::type;
        return std::get<ReacID>(reacM2C_.at(reac));
    }

    template <typename MCurrID>
    typename ModID2ContID<MCurrID>::curr_type getCurrIdx(MCurrID curr) const {
        using CurrID = typename ModID2ContID<MCurrID>::curr_type;
        return std::get<CurrID>(currM2C_.at(curr));
    }

    template <typename rdefT, typename MReacID>
    rdefT& getReacdef(MReacID reac) {
        return *reacdefs<rdefT>().at(getReacIdx(reac).get());
    }

    const Compdef& getInnerComp() const noexcept;

    inline const std::optional<model::compartment_id>& getOuterCompId() const noexcept {
        return outer_compartment_id_;
    }

    model::species_id getSpecModelIdx(container::species_id species) const;

    inline osh::LO getNSpecs() const noexcept {
        return static_cast<osh::LO>(specC2M_.size());
    }

    inline bool getSpecClamped(container::species_id spec) const noexcept {
        return clamped[spec.get()];
    }

    inline void setSpecClamped(container::species_id spec, bool _clampled) noexcept {
        clamped[spec.get()] = _clampled;
    }

    const SReacdef& getReac(container::surface_reaction_id reaction_id) const;

    /**
     * \return the reaction definitions
     */
    template <typename rdefT>
    inline std::vector<std::unique_ptr<rdefT>>& reacdefs() noexcept {
        if constexpr (std::is_same_v<rdefT, SReacdef>) {
            return reacdefPtrs_;
        } else if constexpr (std::is_same_v<rdefT, ComplexSReacdef>) {
            return complexReacdefPtrs_;
        } else if constexpr (std::is_same_v<rdefT, VDepComplexSReacdef>) {
            return vdepComplexReacdefPtrs_;
        } else if constexpr (std::is_same_v<rdefT, VDepSReacdef>) {
            return vdepSReacPtrs_;
        } else if constexpr (std::is_same_v<rdefT, GHKSReacdef>) {
            return ghkSReacPtrs_;
        } else if constexpr (std::is_same_v<rdefT, ComplexGHKSReacdef>) {
            return complexGhkSReacPtrs_;
        } else {
            static_assert(util::always_false_v<rdefT>, "Unmanaged reaction type");
        }
    }

    template <typename CurrT>
    inline const std::vector<std::unique_ptr<CurrT>>& currents() const noexcept {
        if constexpr (std::is_same_v<CurrT, OhmicCurrdef>) {
            return ohmicCurrPtrs;
        } else if constexpr (std::is_same_v<CurrT, ComplexOhmicCurrdef>) {
            return complexOhmicCurrPtrs;
        } else if constexpr (std::is_same_v<CurrT, GHKCurrdef>) {
            return ghkCurrPtrs;
        } else if constexpr (std::is_same_v<CurrT, ComplexGHKCurrdef>) {
            return complexGhkCurrPtrs;
        } else {
            static_assert(util::always_false_v<CurrT>, "Unmanaged current type");
        }
    }

    template <typename CurrT>
    inline std::vector<std::unique_ptr<CurrT>>& currents() noexcept {
        return const_cast<std::vector<std::unique_ptr<CurrT>>&>(
            const_cast<const Patchdef*>(this)->currents<CurrT>());
    }


    inline const Statedef& statedef() const noexcept {
        return pStatedef_;
    }

    inline const std::set<container::species_id>& getAllSpeciesDiffused() const noexcept {
        return species_diffused_;
    }

    inline bool isDiffused(const container::species_id& species) const {
        return species_diffused_.find(species) != species_diffused_.end();
    }

    inline osh::I64 getNReacs() const;

    inline osh::I64 getNKProcs() const;

    container::species_id getSpecContainerIdx(const steps::model::Spec& spec) const;
    container::species_id getSpecContainerIdx(model::species_name species) const;
    container::species_id getSpecContainerIdx(model::species_id species) const;

    //-------------------------------------------------------

    void addGHKReacs(const steps::model::GHKcurr& curr) {
        addReac(curr, true);
        addReac(curr, false);
    }

    void addGHKReacs(const steps::model::ComplexGHKcurr& curr) {
        addReac(curr, true);
        addReac(curr, false);
    }

    //-------------------------------------------------------

    template <typename CurrT, class... ArgsT>
    void addCurrent(const CurrT& curr, ArgsT... args) {
        using CurrdefT = typename Model2Def<CurrT>::curr_def_type;
        using ModelCurrID = typename Model2Def<CurrT>::model_id_type;
        using ContCurrID = typename ModID2ContID<ModelCurrID>::curr_type;

        ModelCurrID current(curr.getID());
        auto& currdefPtrs = currents<CurrdefT>();
        ContCurrID curr_container_idx(static_cast<osh::I64>(currdefPtrs.size()));
        currdefPtrs.emplace_back(std::make_unique<CurrdefT>(*this, curr, args...));
        currM2C_.emplace(current, curr_container_idx);
    }

  private:
    container::species_id addSpec(model::species_id species);

    template <typename SReacT, class... ArgsT>
    void addReac(const SReacT& sreac, ArgsT... args) {
        using RdefT = typename Model2Def<SReacT>::def_type;
        using ModelReacID = typename Model2Def<SReacT>::model_id_type;
        using ContReacID = typename ModID2ContID<ModelReacID>::type;

        ModelReacID reaction(sreac.getID());
        container::kproc_id kproc_id(nKProcs_);
        auto& rdefPtrs = reacdefs<RdefT>();
        ContReacID reac_container_idx(static_cast<osh::I64>(rdefPtrs.size()));
        rdefPtrs.emplace_back(
            std::make_unique<RdefT>(*this, kproc_id, reac_container_idx, sreac, args...));
        reacM2C_.emplace(reaction, reac_container_idx);
        nKProcs_++;
    }

    const DistPatch& pPatch_;
    // compartment KProc order: Reac then Diff
    const Statedef& pStatedef_;
    model::patch_id model_patch_;
    model::compartment_id inner_compartment_id_;
    std::optional<model::compartment_id> outer_compartment_id_;
    container::patch_id container_patch_id_;
    std::unordered_map<model::species_id, container::species_id> specM2C_;
    std::map<AllModSReacID, AllContSReacID> reacM2C_;
    std::map<AllModCurrentID, AllContCurrentID> currM2C_;
    std::vector<model::species_id> specC2M_;
    osh::I64 nKProcs_;

    std::vector<std::unique_ptr<SReacdef>> reacdefPtrs_;
    std::vector<std::unique_ptr<ComplexSReacdef>> complexReacdefPtrs_;
    std::vector<std::unique_ptr<VDepComplexSReacdef>> vdepComplexReacdefPtrs_;
    std::vector<std::unique_ptr<VDepSReacdef>> vdepSReacPtrs_;
    std::vector<std::unique_ptr<GHKSReacdef>> ghkSReacPtrs_;
    std::vector<std::unique_ptr<ComplexGHKSReacdef>> complexGhkSReacPtrs_;

    std::vector<std::unique_ptr<OhmicCurrdef>> ohmicCurrPtrs;
    std::vector<std::unique_ptr<ComplexOhmicCurrdef>> complexOhmicCurrPtrs;
    std::vector<std::unique_ptr<GHKCurrdef>> ghkCurrPtrs;
    std::vector<std::unique_ptr<ComplexGHKCurrdef>> complexGhkCurrPtrs;

    // This is in preparation for SReac
    std::set<container::species_id> species_diffused_;
    std::vector<bool> clamped;
};

}  // namespace steps::dist
