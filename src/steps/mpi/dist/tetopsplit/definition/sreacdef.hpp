
#pragma once

#include <algorithm>
#include <map>
#include <tuple>
#include <vector>

#include "fwd.hpp"
#include "util/vocabulary.hpp"

namespace steps {
namespace dist {

/**
 * \brief Container for a generic surface reactions.
 * Generally, two types of surface reactions exist.
 *   - One for which the \a reaction \a constant does not change throughout the
 *   simulation. In this case the PropensityType argument of this class is a
 * scalar.
 *   - The other one for which the \a reaction \a constant can vary in
 * accordance with voltage on the surface -- so called \a VDepSReac. In that
 * case, the PropensityType argument of this class is a function which takes a
 * voltage as input and returns a reaction constant.
 *
 * \tparam PropensityInfo The type of propensity information of, defining also,
 * the reaction.
 */
template <typename PropensityInfo> class SReacdefBase {
public:
  enum class PoolChangeType { LHS, UPD };
  enum class SpecieClassifier { Reactant, Product };
  enum class SpecieLocation { Patch, InnerCompartment, OuterCompartment };

  using SurfaceReactionComponent =
      std::tuple<container::species_id, SpecieClassifier, SpecieLocation>;
  using SurfaceReactionComponents = std::vector<SurfaceReactionComponent>;

  SReacdefBase(const Statedef &state_def, const container::patch_id &patch_id,
               container::kproc_id kproc_id,
               container::surface_reaction_id surface_reaction_id,
               const SurfaceReactionComponents &reaction_definition,
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

  inline container::surface_reaction_id surfaceReactionID() const noexcept {
    return surface_reaction_;
  }

  inline const PropensityInfo &getInfo() const noexcept { return info_; }

  inline osh::I64 getOrder() const noexcept { return order_; }

  const Patchdef &patchdef() const noexcept;

  template <PoolChangeType PoolType, SpecieLocation Location>
  inline const std::map<container::species_id, osh::LO> &
  getStoichiometry() const noexcept;

private:
  const Statedef &state_def_;
  const container::patch_id patch_id_;

  const container::kproc_id kproc_id_;
  const container::surface_reaction_id surface_reaction_;
  const SurfaceReactionComponents &reaction_components_;

  std::map<container::species_id, osh::LO> patch_lhs_;
  std::map<container::species_id, osh::LO> patch_update_;
  std::map<container::species_id, osh::LO> inner_lhs_;
  std::map<container::species_id, osh::LO> inner_update_;
  std::map<container::species_id, osh::LO> outer_lhs_;
  std::map<container::species_id, osh::LO> outer_update_;

  osh::I64 order_;
  const PropensityInfo info_;
  bool is_inner_, is_surf_surf_;
};

//-------------------------------------------------------

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::LHS,
                              GHKSReacdef::SpecieLocation::Patch>()
    const noexcept {
  return patch_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS,
                           SReacdef::SpecieLocation::Patch>() const noexcept {
  return patch_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<VDepSReacdef::PoolChangeType::LHS,
                               VDepSReacdef::SpecieLocation::Patch>() const
    noexcept {
  return patch_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS,
                           SReacdef::SpecieLocation::InnerCompartment>() const
    noexcept {
  return inner_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<VDepSReacdef::PoolChangeType::LHS,
                               VDepSReacdef::SpecieLocation::InnerCompartment>()
    const noexcept {
  return inner_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::LHS,
                              GHKSReacdef::SpecieLocation::InnerCompartment>()
    const noexcept {
  return inner_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::LHS,
                           SReacdef::SpecieLocation::OuterCompartment>() const
    noexcept {
  return outer_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<VDepSReacdef::PoolChangeType::LHS,
                               VDepSReacdef::SpecieLocation::OuterCompartment>()
    const noexcept {
  return outer_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::LHS,
                              GHKSReacdef::SpecieLocation::OuterCompartment>()
    const noexcept {
  return outer_lhs_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD,
                           SReacdef::SpecieLocation ::Patch>() const noexcept {
  return patch_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<VDepSReacdef::PoolChangeType::UPD,
                               VDepSReacdef::SpecieLocation ::Patch>() const
    noexcept {
  return patch_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::UPD,
                              GHKSReacdef::SpecieLocation ::Patch>()
    const noexcept {
  return patch_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD,
                           SReacdef::SpecieLocation ::InnerCompartment>() const
    noexcept {
  return inner_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<
    VDepSReacdef::PoolChangeType::UPD,
    VDepSReacdef::SpecieLocation ::InnerCompartment>() const noexcept {
  return inner_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::UPD,
                              GHKSReacdef::SpecieLocation ::InnerCompartment>()
    const noexcept {
  return inner_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
SReacdef::getStoichiometry<SReacdef::PoolChangeType::UPD,
                           SReacdef::SpecieLocation ::OuterCompartment>() const
    noexcept {
  return outer_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
VDepSReacdef::getStoichiometry<
    VDepSReacdef::PoolChangeType::UPD,
    VDepSReacdef::SpecieLocation ::OuterCompartment>() const noexcept {
  return outer_update_;
}

template <>
template <>
inline const std::map<container::species_id, osh::LO> &
GHKSReacdef::getStoichiometry<GHKSReacdef::PoolChangeType::UPD,
                              GHKSReacdef::SpecieLocation ::OuterCompartment>()
    const noexcept {
  return outer_update_;
}

//-------------------------------------------------------

extern template class SReacdefBase<GHKInfo>;
extern template class SReacdefBase<SReacInfo>;
extern template class SReacdefBase<VDepInfo>;

}  // namespace dist
} // namespace steps
