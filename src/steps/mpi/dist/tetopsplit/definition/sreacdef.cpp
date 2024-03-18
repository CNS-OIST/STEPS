#include "sreacdef.hpp"

#include "patchdef.hpp"
#include "statedef.hpp"

namespace steps::dist {

//-------------------------------------------------------

template <typename PropensityInfo>
SReacdefBase<PropensityInfo>::SReacdefBase(const Statedef& state_def,
                                           const container::patch_id& patch_id,
                                           container::kproc_id kproc_id,
                                           container::surface_reaction_id surface_reaction_id,
                                           const SurfaceReactionComponents& reaction_definition,
                                           PropensityInfo info)
    : state_def_(state_def)
    , patch_id_(patch_id)
    , kproc_id_(kproc_id)
    , surface_reaction_(surface_reaction_id)
    , reaction_components_(reaction_definition)
    , info_(info) {
    order_ = 0;
    for (const auto& t: reaction_components_) {
        const auto& [species_id, speciesClassifier, surfaceReactionSpecieLocation] = t;
        if (speciesClassifier == SpecieClassifier::Reactant) {
            order_++;
        }
        auto populateStoichiometry = [](std::map<container::species_id, osh::LO>& lhs,
                                        std::map<container::species_id, osh::LO>& upd,
                                        SpecieClassifier classifier,
                                        container::species_id species) {
            if (classifier == SpecieClassifier::Reactant) {
                auto itlhs = lhs.find(species);
                if (itlhs == lhs.end()) {
                    lhs[species] = -1;
                } else {
                    itlhs->second--;
                }
                auto itupd = upd.find(species);
                if (itupd == upd.end()) {
                    upd[species] = -1;
                } else {
                    itupd->second--;
                }
            } else {
                assert(classifier == SpecieClassifier::Product);
                auto itrhs = upd.find(species);
                if (itrhs == upd.end()) {
                    upd[species] = 1;
                } else {
                    itrhs->second++;
                }
            }
        };
        if (surfaceReactionSpecieLocation == SpecieLocation::InnerCompartment) {
            populateStoichiometry(inner_lhs_, inner_update_, speciesClassifier, species_id);
        }
        if (surfaceReactionSpecieLocation == SpecieLocation::Patch) {
            populateStoichiometry(patch_lhs_, patch_update_, speciesClassifier, species_id);
        }
        if (surfaceReactionSpecieLocation == SpecieLocation::OuterCompartment) {
            populateStoichiometry(outer_lhs_, outer_update_, speciesClassifier, species_id);
        }
    }

    auto removeZeroElements = [](std::map<container::species_id, osh::LO>& d) {
        auto it = d.begin();
        while (it != d.end()) {
            if (it->second == 0) {
                it = d.erase(it);
            } else {
                it++;
            }
        }
    };
    removeZeroElements(inner_update_);
    removeZeroElements(outer_update_);
    removeZeroElements(patch_update_);

    const bool all_inner = outer_lhs_.empty();
    const bool all_outer = inner_lhs_.empty();
    if (!(all_inner || all_outer) &&
        typeid(PropensityInfo).hash_code() != typeid(GHKInfo).hash_code()) {
        throw std::logic_error(
            "A surface reaction involves volume reactants "
            "all in the inner or all in the outer compartment");
    }
    is_inner_ = all_inner;
    is_surf_surf_ = outer_lhs_.empty() && inner_lhs_.empty();
}

//-------------------------------------------------------

template <typename PropensityInfo>
const Patchdef& SReacdefBase<PropensityInfo>::patchdef() const noexcept {
    return state_def_.getPatchdef(patch_id_);
}

//-------------------------------------------------------

template class SReacdefBase<GHKInfo>;
template class SReacdefBase<SReacInfo>;
template class SReacdefBase<VDepInfo>;

}  // namespace steps::dist
