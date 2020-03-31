
#include "opsplit/patchdef.hpp"
#include "opsplit/sreacdef.hpp"
#include "opsplit/statedef.hpp"

namespace zee {

//-------------------------------------------------------

Patchdef::Patchdef(const Statedef& statedef,
                   model::patch_id t_model_patch,
                   model::compartment_id inner_compartment_id,
                   const boost::optional<model::compartment_id>& outer_compartment_id)
    : pStatedef_(statedef)
    , model_patch_(t_model_patch)
    , inner_compartment_id_(inner_compartment_id)
    , outer_compartment_id_(outer_compartment_id)
    , nKProcs_(0){};

//-------------------------------------------------------

container::surface_reaction_id Patchdef::addSurfaceReac(
    const std::vector<container::specie_id>& reactants_i,
    const std::vector<container::specie_id>& reactants_s,
    const std::vector<container::specie_id>& reactants_o,
    const std::vector<container::specie_id>& product_i,
    const std::vector<container::specie_id>& product_s,
    const std::vector<container::specie_id>& product_o,
    PetscScalar kcst) {
    SReacdef::SurfaceReactionComponents reaction_components;
    reaction_components.reserve(reactants_i.size() + reactants_s.size() + reactants_o.size() +
                                product_i.size() + product_s.size() + product_o.size());

    // prepare arguments for building SReacdef
    auto accum_reaction_comps = [&reaction_components,
                                 this](const auto loc, const auto classifier, const auto& vec) {
        return std::transform(vec.begin(),
                              vec.end(),
                              std::back_inserter(reaction_components),
                              [&loc, &classifier, this](const auto& v) {
                                  if (loc == SReacdef::SpecieLocation::Patch) {
                                      getSpecModelIdx(v);  // Test whether specie is registered
                                  }
                                  return std::make_tuple(v, classifier, loc);
                              });
    };
    accum_reaction_comps(SReacdef::SpecieLocation::InnerCompartment,
                         SReacdef::SpecieClassifier::Reactant,
                         reactants_i);
    accum_reaction_comps(SReacdef::SpecieLocation::Patch,
                         SReacdef::SpecieClassifier::Reactant,
                         reactants_s);
    accum_reaction_comps(SReacdef::SpecieLocation::OuterCompartment,
                         SReacdef::SpecieClassifier::Reactant,
                         reactants_o);
    accum_reaction_comps(SReacdef::SpecieLocation::InnerCompartment,
                         SReacdef::SpecieClassifier::Product,
                         product_i);
    accum_reaction_comps(SReacdef::SpecieLocation::Patch,
                         SReacdef::SpecieClassifier::Product,
                         product_s);
    accum_reaction_comps(SReacdef::SpecieLocation::OuterCompartment,
                         SReacdef::SpecieClassifier::Product,
                         product_o);

    container::surface_reaction_id reac_id(static_cast<PetscInt>(reacdefPtrs_.size()));
    reacdefPtrs_.push_back(std::make_unique<SReacdef>(pStatedef_,
                                                      container_patch_,
                                                      container::kproc_id(nKProcs_),
                                                      reac_id,
                                                      reaction_components,
                                                      kcst));
    // nkprocs_ records the number of kprocs in the patch
    nKProcs_++;
    return reac_id;
}

//-------------------------------------------------------

model::specie_id Patchdef::getSpecModelIdx(container::specie_id specie) const {
    if (!(specie < static_cast<container::specie_id::value_type>(specC2M_.size()))) {
        throw std::invalid_argument(std::string("Unregistered needed specie in patch ") +
                                    model_patch_);
    }
    return specC2M_[static_cast<size_t>(specie.get())];
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

container::specie_id Patchdef::addSpec(model::specie_id specie) {
    auto specieIt = specM2C_.find(specie);
    if (specieIt != specM2C_.end()) {
        return specieIt->second;
    }
    const container::specie_id spec_container_idx(
        static_cast<container::specie_id::value_type>(specC2M_.size()));
    specM2C_[specie] = spec_container_idx;
    specC2M_.push_back(specie);
    return spec_container_idx;
}

//-------------------------------------------------------

inline PetscInt Patchdef::getNReacs() const {
    return static_cast<PetscInt>(reacdefPtrs_.size());
}

//-------------------------------------------------------

inline KProcType Patchdef::getKProcType(container::kproc_id kproc) const {
    if (kproc.get() >= 0 && kproc < getNReacs()) {
        return KProcType::SReac;
    } else {
        throw std::out_of_range("KProc local index error.");
    }
}

//-------------------------------------------------------

inline PetscInt Patchdef::getNKProcs() const {
    return nKProcs_;
}

//-------------------------------------------------------
}  // namespace zee
