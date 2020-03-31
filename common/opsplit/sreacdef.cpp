#include "opsplit/sreacdef.hpp"
#include "opsplit/patchdef.hpp"
#include "opsplit/statedef.hpp"

namespace zee {

//-------------------------------------------------------

SReacdef::SReacdef(const Statedef& state_def,
                   const container::patch_id patch_id,
                   const container::kproc_id kproc_id,
                   const container::surface_reaction_id surface_reaction_id,
                   const SurfaceReactionComponents& reaction_definition,
                   const PetscScalar kcst)
    : state_def_(state_def)
    , patch_id_(patch_id)
    , kproc_id_(kproc_id)
    , surface_reaction_(surface_reaction_id)
    , reaction_components_(reaction_definition)
    , kcst_(kcst) {
    container::specie_id specie_id;
    SpecieClassifier specieClassifier;
    SpecieLocation surfaceReactionSpecieLocation;
    order_ = 0;
    for (const auto& t: reaction_components_) {
        std::tie(specie_id, specieClassifier, surfaceReactionSpecieLocation) = t;
        if (specieClassifier == SpecieClassifier::Reactant) {
            order_++;
        }
        auto populateStoichiometry = [](std::map<container::specie_id, osh::LO>& lhs,
                                        std::map<container::specie_id, osh::LO>& upd,
                                        SpecieClassifier classifier,
                                        container::specie_id specie) {
            if (classifier == SpecieClassifier::Reactant) {
                auto itlhs = lhs.find(specie);
                if (itlhs == lhs.end()) {
                    lhs[specie] = -1;
                } else {
                    itlhs->second--;
                }
                auto itupd = upd.find(specie);
                if (itupd == upd.end()) {
                    upd[specie] = -1;
                } else {
                    (*itupd).second--;
                }
            } else {
                assert(classifier == SpecieClassifier::Product);
                auto itrhs = upd.find(specie);
                if (itrhs == upd.end()) {
                    upd[specie] = 1;
                } else {
                    (*itrhs).second++;
                }
            }
        };
        if (surfaceReactionSpecieLocation == SpecieLocation::InnerCompartment) {
            populateStoichiometry(inner_lhs_, inner_update_, specieClassifier, specie_id);
        }
        if (surfaceReactionSpecieLocation == SpecieLocation::Patch) {
            populateStoichiometry(patch_lhs_, patch_update_, specieClassifier, specie_id);
        }
        if (surfaceReactionSpecieLocation == SpecieLocation::OuterCompartment) {
            populateStoichiometry(outer_lhs_, outer_update_, specieClassifier, specie_id);
        }
    }

    bool all_inner = outer_lhs_.size() == 0;
    bool all_outer = inner_lhs_.size() == 0;
    if (!(all_inner || all_outer)) {
        throw std::logic_error(
            "A surface reaction involves reactants volume reactants"
            "all in the inner or all in the outer compartment");
    }
    is_inner_ = all_inner;
    is_surf_surf_ = (outer_lhs_.size() == 0) && (inner_lhs_.size() == 0);
}

//-------------------------------------------------------

const Patchdef& SReacdef::patchdef() const {
    return state_def_.getPatchdef(patch_id_);
}

//-------------------------------------------------------

};  // namespace zee
