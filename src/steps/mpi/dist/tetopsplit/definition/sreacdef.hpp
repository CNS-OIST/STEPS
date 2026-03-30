
#pragma once

#include <algorithm>
#include <map>
#include <tuple>
#include <type_traits>
#include <vector>

#include "compdef.hpp"
#include "complexeventsdef.hpp"
#include "fwd.hpp"
#include "math/constants.hpp"
#include "model/complexsreac.hpp"
#include "model/ghkcurr.hpp"
#include "model/vdepsreac.hpp"
#include "patchdef.hpp"
#include "statedef.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist {

enum class PoolChangeType { LHS, UPD };
enum class SpecieClassifier { Reactant, Product };

using SurfaceReactionComponent =
    std::tuple<container::species_id, SpecieClassifier, steps::model::Location>;
using SurfaceReactionComponents = std::vector<SurfaceReactionComponent>;

/**
 * getPropInfo functions take a steps::model reaction or current object as parameter and return the
 * corresponding PropensityInfo.
 */
template <typename VDepSReacT,
          typename = std::enable_if_t<std::is_same_v<VDepSReacT, steps::model::VDepSReac> ||
                                      std::is_same_v<VDepSReacT, steps::model::VDepComplexSReac>>>
VDepInfo getPropInfo(const VDepSReacT& reac);
template <typename CurrType,
          typename = std::enable_if_t<std::is_same_v<CurrType, steps::model::GHKcurr> ||
                                      std::is_same_v<CurrType, steps::model::ComplexGHKcurr>>>
GHKInfo getPropInfo(const CurrType& curr, bool in2out = true);
template <typename SReacT,
          typename = std::enable_if_t<std::is_same_v<SReacT, steps::model::SReac> ||
                                      std::is_same_v<SReacT, steps::model::ComplexSReac>>>
SReacInfo getPropInfo(const SReacT& reac);

/**
 * getReactionComponents functions take a steps::model reaction or current object as parameter and
 * return a SurfaceReactionComponents that contains the species involved in the reaction.
 */
SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef,
                                                const steps::model::GHKcurr& curr,
                                                bool in2out);
SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef,
                                                const steps::model::ComplexGHKcurr& curr,
                                                bool in2out);
template <typename SReacT>
SurfaceReactionComponents getReactionComponents(const Patchdef& patchdef, const SReacT& reac);

// Mapping between steps::model reaction or current type and PropensityInfo type
template <typename SReacT>
using Model2PropInfo = decltype(getPropInfo(std::declval<const SReacT&>()));

/**
 * \brief Container for a generic surface reactions.
 * Generally, two types of surface reactions exist.
 *   - One for which the \a reaction \a constant does not change throughout the
 *   simulation. In this case the PropensityInfo argument of this class is a
 * scalar.
 *   - The other one for which the \a reaction \a constant can vary in
 * accordance with voltage on the surface -- so called \a VDepSReac. In that
 * case, the PropensityInfo argument of this class is a function which takes a
 * voltage as input and returns a reaction constant.
 *
 * \tparam PropensityInfo The type of propensity information of, defining also,
 * the reaction.
 */
template <typename SReacT>
class SReacdefBase {
  public:
    using PropensityInfo = Model2PropInfo<SReacT>;
    using ContReacID = typename ModID2ContID<typename Model2Def<SReacT>::model_id_type>::type;

    SReacdefBase(const Patchdef& patchdef,
                 container::kproc_id kproc,
                 ContReacID surface_reaction_id,
                 const SurfaceReactionComponents& reaction_definition,
                 PropensityInfo info);

    bool isInnerCompartmentReaction() const noexcept {
        // reactant volume species need to be either in the inner compartment
        // or the outer compartment as per STEPS
        return is_inner_;
    }

    inline bool isSurfaceSurfaceReaction() const noexcept {
        // check whether all reactants are on the surface
        return is_surf_surf_;
    }

    inline container::kproc_id getKProcContainerIdx() const noexcept {
        return kproc_id_;
    }

    inline ContReacID surfaceReactionID() const noexcept {
        return surface_reaction_;
    }

    inline PropensityInfo& getInfo() noexcept {
        return info_;
    }

    inline osh::I64 getOrder() const noexcept {
        return order_;
    }

    const Patchdef& patchdef() const noexcept {
        return patchdef_;
    }

    inline osh::I64 getExtent() const noexcept {
        return extent;
    }

    inline void add_extent(int nb) noexcept {
        extent += nb;
    }

    inline const std::map<container::species_id, osh::LO>& getLHS(
        steps::model::Location loc) const {
        return lhs_.at(loc);
    }
    inline const std::map<container::species_id, osh::LO>& getUPD(
        steps::model::Location loc) const {
        return upd_.at(loc);
    }

    template <PoolChangeType tpe>
    inline const std::map<container::species_id, osh::LO>& getStoichiometry(
        steps::model::Location loc) const {
        if constexpr (tpe == PoolChangeType::LHS) {
            return lhs_.at(loc);
        } else if constexpr (tpe == PoolChangeType::UPD) {
            return upd_.at(loc);
        } else {
            static_assert(tpe != tpe, "Unexpected PoolChange value.");
        }
    }

    void reset() {
        extent = 0;
    }

  protected:
    const Patchdef& patchdef_;

    const container::kproc_id kproc_id_;
    const ContReacID surface_reaction_;

    std::map<steps::model::Location, std::map<container::species_id, osh::LO>> lhs_;
    std::map<steps::model::Location, std::map<container::species_id, osh::LO>> upd_;

    osh::I64 order_{0};
    osh::I64 extent{0};
    PropensityInfo info_;
    bool is_inner_, is_surf_surf_;
};

/**
 * \brief Generic class for surface reactions that are declared as reactions in steps::model
 *
 * \tparam SReacT The steps::model type of the surface reaction
 */
template <typename SReacT>
class ModelSReacdef: public SReacdefBase<SReacT> {
  public:
    using ContReacID = typename SReacdefBase<SReacT>::ContReacID;

    ModelSReacdef(const Patchdef& patchdef,
                  container::kproc_id kproc,
                  ContReacID surface_reaction_id,
                  const SReacT& reac);

    /// Reset the def-object values to model defaults
    void reset();

  protected:
    const SReacT& reac_;
};

/**
 * \brief Base class for GHK current surface reaction
 */
template <typename GHKcurrT>
class GHKSReacdefBase: public SReacdefBase<GHKcurrT> {
  public:
    using ContReacID = typename SReacdefBase<GHKcurrT>::ContReacID;

    GHKSReacdefBase(const Patchdef& patchdef,
                    container::kproc_id kproc,
                    ContReacID surface_reaction_id,
                    const GHKcurrT& curr,
                    bool in2out);

    /// Reset the def-object values to model defaults
    void reset();

  protected:
    const GHKcurrT& curr_;
    bool in2out_;
};

/**
 * \brief Surface reaction for GHK current with species-like channel state
 */
class GHKSReacdef: public GHKSReacdefBase<steps::model::GHKcurr> {
  public:
    using ContReacID = typename GHKSReacdefBase<steps::model::GHKcurr>::ContReacID;

    GHKSReacdef(const Patchdef& patchdef,
                container::kproc_id kproc,
                ContReacID surface_reaction_id,
                const steps::model::GHKcurr& curr,
                bool in2out);
};

/**
 * Surface reaction for GHK current with complex channel state
 *
 * Only one complex update events can be involved since the current cannot lead to creation or
 * deletion of channels.
 */
class ComplexGHKSReacdef: public GHKSReacdefBase<steps::model::ComplexGHKcurr> {
  public:
    using ContReacID = typename GHKSReacdefBase<steps::model::ComplexGHKcurr>::ContReacID;

    ComplexGHKSReacdef(const Patchdef& patchdef,
                       container::kproc_id kproc,
                       ContReacID surface_reaction_id,
                       const steps::model::ComplexGHKcurr& curr,
                       bool in2out);

    const std::set<model::complex_substate_id>& complexDEPSET() const noexcept {
        return pComplex_DEPSET;
    }

    const std::shared_ptr<ComplexUpdateEventdef>& updEvent() const noexcept {
        return pComplexSurfUPDEv;
    }

    model::complex_id complexID() const {
        return pComplexId;
    }

  private:
    model::complex_id pComplexId;
    std::shared_ptr<ComplexUpdateEventdef> pComplexSurfUPDEv;

    std::set<model::complex_substate_id> pComplex_DEPSET;
};

/**
 * \brief Generic class for complex surface reactions that are declared as reactions in steps::model
 *
 * \tparam CSReacT The steps::model type of the complex surface reaction
 */
template <typename CSReacT>
class ModelComplexSReacdef: public ModelSReacdef<CSReacT> {
  public:
    using ContReacID = typename ModelSReacdef<CSReacT>::ContReacID;

    ModelComplexSReacdef(const Patchdef& patchdef,
                         container::kproc_id kproc,
                         ContReacID surface_reaction_id,
                         const CSReacT& reac);

    const std::map<steps::model::Location,
                   std::map<model::complex_id, std::set<model::complex_substate_id>>>&
    complexDEPMAP() const noexcept {
        return pComplex_DEPMAP;
    }

    const std::map<steps::model::Location,
                   std::map<model::complex_id, std::set<model::complex_substate_id>>>&
    complexUPDMAP() const noexcept {
        return pComplex_UPDMAP;
    }

    const std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexUpdateEventdef>>>&
    updEvents() const noexcept {
        return pComplexUPDEvs;
    }

    const std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexDeleteEventdef>>>&
    delEvents() const noexcept {
        return pComplexDELEvs;
    }

    const std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexCreateEventdef>>>&
    creEvents() const noexcept {
        return pComplexCREEvs;
    }

    std::vector<std::shared_ptr<ComplexLHSEventdef>> lhsEvents(steps::model::Location loc) const {
        std::vector<std::shared_ptr<ComplexLHSEventdef>> ret;
        ret.reserve(pComplexUPDEvs.at(loc).size() + pComplexDELEvs.at(loc).size());
        auto& updevs = pComplexUPDEvs.at(loc);
        auto& delevs = pComplexDELEvs.at(loc);
        ret.insert(ret.end(), updevs.begin(), updevs.end());
        ret.insert(ret.end(), delevs.begin(), delevs.end());
        return ret;
    }

  private:
    std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexUpdateEventdef>>>
        pComplexUPDEvs;
    std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexDeleteEventdef>>>
        pComplexDELEvs;
    std::map<steps::model::Location, std::vector<std::shared_ptr<ComplexCreateEventdef>>>
        pComplexCREEvs;

    // location -> {cmplxIdx -> {sub unit states ind}}
    std::map<steps::model::Location,
             std::map<model::complex_id, std::set<model::complex_substate_id>>>
        pComplex_DEPMAP;
    std::map<steps::model::Location,
             std::map<model::complex_id, std::set<model::complex_substate_id>>>
        pComplex_UPDMAP;
};

//-------------------------------------------------------

}  // namespace steps::dist
