
#include "patchdef.hpp"

#include "sreacdef.hpp"
#include "statedef.hpp"

namespace steps::dist {

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<SReacdef>>& Patchdef::getContainer<SReacInfo>() {
    return reacdefPtrs_;
}

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<VDepSReacdef>>& Patchdef::getContainer<VDepInfo>() {
    return vdepSReacPtrs_;
}

//-------------------------------------------------------

template <>
std::vector<std::unique_ptr<GHKSReacdef>>& Patchdef::getContainer<GHKInfo>() {
    return ghkSReacPtrs_;
}

//-------------------------------------------------------

container::species_id Patchdef::getSpecPatchIdx(model::species_id species) const {
    const auto& it = specM2C_.find(species);
    if (it != specM2C_.end()) {
        return it->second;
    }
    throw std::logic_error("Unregistered species id " + pStatedef_.getSpecID(species));
}

//-------------------------------------------------------

Patchdef::Patchdef(const Statedef& statedef,
                   model::patch_id t_model_patch,
                   container::patch_id container_patch_id,
                   model::compartment_id inner_compartment_id,
                   const std::optional<model::compartment_id>& outer_compartment_id)
    : pStatedef_(statedef)
    , model_patch_(t_model_patch)
    , inner_compartment_id_(inner_compartment_id)
    , outer_compartment_id_(outer_compartment_id)
    , container_patch_id_(container_patch_id)
    , nKProcs_(0) {}

//-------------------------------------------------------

template <typename PropensityType>
container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id>& reactants_i,
    const std::vector<container::species_id>& reactants_s,
    const std::vector<container::species_id>& reactants_o,
    const std::vector<container::species_id>& products_i,
    const std::vector<container::species_id>& products_s,
    const std::vector<container::species_id>& products_o,
    PropensityType kcst) {
    typename SReacdefBase<PropensityType>::SurfaceReactionComponents reaction_components;
    reaction_components.reserve(reactants_i.size() + reactants_s.size() + reactants_o.size() +
                                products_i.size() + products_s.size() + products_o.size());

    // prepare arguments for building SReacdef
    auto accum_reaction_comps = [&reaction_components,
                                 this](const auto loc, const auto classifier, const auto& vec) {
        return std::transform(vec.begin(),
                              vec.end(),
                              std::back_inserter(reaction_components),
                              [&loc, &classifier, this](const auto& v) {
                                  if (loc == SReacdefBase<PropensityType>::SpecieLocation::Patch) {
                                      getSpecModelIdx(v);  // Test whether species is registered
                                  }
                                  return std::make_tuple(v, classifier, loc);
                              });
    };
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::InnerCompartment,
                         SReacdefBase<PropensityType>::SpecieClassifier::Reactant,
                         reactants_i);
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::Patch,
                         SReacdefBase<PropensityType>::SpecieClassifier::Reactant,
                         reactants_s);
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::OuterCompartment,
                         SReacdefBase<PropensityType>::SpecieClassifier::Reactant,
                         reactants_o);
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::InnerCompartment,
                         SReacdefBase<PropensityType>::SpecieClassifier::Product,
                         products_i);
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::Patch,
                         SReacdefBase<PropensityType>::SpecieClassifier::Product,
                         products_s);
    accum_reaction_comps(SReacdefBase<PropensityType>::SpecieLocation::OuterCompartment,
                         SReacdefBase<PropensityType>::SpecieClassifier::Product,
                         products_o);

    container::surface_reaction_id reac_id(
        static_cast<osh::I64>(getContainer<PropensityType>().size()));

    getContainer<PropensityType>().push_back(
        std::make_unique<SReacdefBase<PropensityType>>(pStatedef_,
                                                       container_patch_id_,
                                                       container::kproc_id(nKProcs_),
                                                       reac_id,
                                                       reaction_components,
                                                       kcst));
    // nkprocs_ records the number of kprocs in the patch
    nKProcs_++;
    return reac_id;
}

//-------------------------------------------------------

model::species_id Patchdef::getSpecModelIdx(container::species_id species) const {
    if (!(species < static_cast<container::species_id::value_type>(specC2M_.size()))) {
        throw std::invalid_argument(std::string("Unregistered needed species in patch ") +
                                    model_patch_);
    }
    return specC2M_[static_cast<size_t>(species.get())];
}

//-------------------------------------------------------

inline const Compdef& Patchdef::getInnerComp() const noexcept {
    return pStatedef_.getCompdef(inner_compartment_id_);
}

//-------------------------------------------------------

const SReacdef& Patchdef::getReac(container::surface_reaction_id reaction_id) const {
    return *reacdefPtrs_[static_cast<size_t>(reaction_id.get())];
}

//-------------------------------------------------------

container::species_id Patchdef::addSpec(model::species_id species) {
    auto speciesIt = specM2C_.find(species);
    if (speciesIt != specM2C_.end()) {
        return speciesIt->second;
    }
    const container::species_id spec_container_idx(
        static_cast<container::species_id::value_type>(specC2M_.size()));
    specM2C_[species] = spec_container_idx;
    specC2M_.push_back(species);
    return spec_container_idx;
}

//-------------------------------------------------------

inline osh::I64 Patchdef::getNReacs() const {
    return static_cast<osh::I64>(reacdefPtrs_.size());
}

//-------------------------------------------------------

inline kproc::KProcType Patchdef::getKProcType(container::kproc_id kproc) const {
    if (kproc.get() >= 0 && kproc < getNReacs()) {
        return kproc::KProcType::SReac;
    } else {
        throw std::out_of_range("KProc local index error.");
    }
}

//-------------------------------------------------------

inline osh::I64 Patchdef::getNKProcs() const {
    return nKProcs_;
}

//-------------------------------------------------------

template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id>& reactants_i,
    const std::vector<container::species_id>& reactants_s,
    const std::vector<container::species_id>& reactants_o,
    const std::vector<container::species_id>& product_i,
    const std::vector<container::species_id>& product_s,
    const std::vector<container::species_id>& product_o,
    SReacInfo kcst);

//-------------------------------------------------------

template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id>& reactants_i,
    const std::vector<container::species_id>& reactants_s,
    const std::vector<container::species_id>& reactants_o,
    const std::vector<container::species_id>& product_i,
    const std::vector<container::species_id>& product_s,
    const std::vector<container::species_id>& product_o,
    VDepInfo kcst);
//-------------------------------------------------------

template container::surface_reaction_id Patchdef::addSurfaceReacImpl(
    const std::vector<container::species_id>& reactants_i,
    const std::vector<container::species_id>& reactants_s,
    const std::vector<container::species_id>& reactants_o,
    const std::vector<container::species_id>& product_i,
    const std::vector<container::species_id>& product_s,
    const std::vector<container::species_id>& product_o,
    GHKInfo kcst);

//-------------------------------------------------------

}  // namespace steps::dist
