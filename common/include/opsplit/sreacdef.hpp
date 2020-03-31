
#pragma once

#include <algorithm>
#include <map>
#include <petscsys.h>
#include <tuple>
#include <vector>

#include "fwd.hpp"
#include "opsplit/vocabulary.hpp"

namespace zee {

class SReacdef {
  public:
    enum class PoolChangeType { LHS, UPD };

    enum class SpecieClassifier { Reactant, Product };
    enum class SpecieLocation { Patch, InnerCompartment, OuterCompartment };

    typedef std::tuple<container::specie_id, SpecieClassifier, SpecieLocation>
        SurfaceReactionComponent;
    typedef std::vector<SurfaceReactionComponent> SurfaceReactionComponents;

    SReacdef(const Statedef& state_def,
             container::patch_id patch_id,
             container::kproc_id kproc_id,
             container::surface_reaction_id surface_reaction_id,
             const SurfaceReactionComponents& reaction_definition,
             PetscScalar kcst);

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

    inline PetscScalar getKcst() const noexcept {
        return kcst_;
    }

    inline PetscInt getOrder() const noexcept {
        return order_;
    }

    const Patchdef& patchdef() const;

    template <PoolChangeType PoolType, SpecieLocation Location>
    inline const std::map<container::specie_id, osh::LO>& getStoichiometry() const;

  private:
    const Statedef& state_def_;
    container::patch_id patch_id_;

    container::kproc_id kproc_id_;
    container::surface_reaction_id surface_reaction_;
    SurfaceReactionComponents reaction_components_;

    std::map<container::specie_id, osh::LO> patch_lhs_;
    std::map<container::specie_id, osh::LO> patch_update_;
    std::map<container::specie_id, osh::LO> inner_lhs_;
    std::map<container::specie_id, osh::LO> inner_update_;
    std::map<container::specie_id, osh::LO> outer_lhs_;
    std::map<container::specie_id, osh::LO> outer_update_;

    PetscInt order_;
    PetscScalar kcst_;
    bool is_inner_, is_surf_surf_;
};

//-------------------------------------------------------

template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS, SReacdef::SpecieLocation::Patch>() const {
    return patch_lhs_;
}
template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS,
                           SReacdef::SpecieLocation::InnerCompartment>() const {
    return inner_lhs_;
}
template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS,
                           SReacdef::SpecieLocation::OuterCompartment>() const {
    return outer_lhs_;
}

template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD, SReacdef::SpecieLocation ::Patch>()
    const {
    return patch_update_;
}
template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD,
                           SReacdef::SpecieLocation ::InnerCompartment>() const {
    return inner_update_;
}
template <>
inline const std::map<container::specie_id, osh::LO>&
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD,
                           SReacdef::SpecieLocation ::OuterCompartment>() const {
    return outer_update_;
}

//-------------------------------------------------------
};  // namespace zee
