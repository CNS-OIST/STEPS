#pragma once

#include <random>

#include <boost/optional/optional.hpp>

#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/propensities.hpp"
#include "mpi/dist/tetopsplit/kproc/reactions.hpp"
#include "rng/rng.hpp"
#include "util/collections.hpp"

namespace steps {
namespace dist {

/**
 * \brief RSSA operator.
 *
 * This class implements the RSSA operator in the OpSplit solution.
 *
 * The core function of this class is run(), which partially corresponds to the
 * SSA component the simulation core loop tetexact::Tetexact::run() /
 * tetopsplit::TetOpSplit::run() in STEPS.
 */
template <typename RNG, typename NumMolecules>
class RSSAOperator {
public:
  template <int Policy>
  using propensities_groups_t =
      kproc::PropensitiesGroup<NumMolecules, Policy>;
  using Event = std::pair<osh::Real, kproc::KProcID>;
  /**
   * Constructor
   * \param mol_state molecular state
   * \param dv molecules transferred to neighbours
   * \param kproc_state kinetic processes state
   */
  RSSAOperator(MolState<NumMolecules>& mol_state,
               kproc::KProcState& kproc_state,
               RNG& t_rng,
               osh::Reals potential_on_vertices);

  /// \return the number of kinetic events that arose on the current process
  inline osh::I64 getExtent() const noexcept { return extent; }

  osh::Real run(osh::Real period, osh::Real state_time);
 
 /* TODO */
 /**
   * \brief reset the required propensity data stored in the operator
   */
  void reset() {
    // make it do nothing now to pass the CI.
  }

 /* TODO */
  /**
   * \brief Set the maximum time to be considered as valid in the SSA.
   * used in GB optimization
   * \param max_time maximum simulation end time for the SSA search system
   */
  void updateMaxTime(const osh::Real /*max_time*/) {
     // make it do nothing now to pass the CI.
  }

private:
  /** Updates bounds for propensity update
   *
   * In the RSSA method propensities are updated only if the molecule count goes out of bound.
   * In that case new bounds are recomputed (using this function). Bounds are there to keep
   * a good Rejection-Acceptance ratio.
   *
   * \tparam Entity: a strong_id type, \c mesh::element_id or \c mesh::triangle_id_t for instance
   * \param num_molecules: number of molecules
   * \param num_molecules_lower_bound: if num_molecules goes below, recompute propensities
   * \param num_molecules_upper_bound: if num_molecules goes above, recompute propensities
   */
  template <typename Entity>
  void updateNumMolsBounds(const EntityMolecules<Entity, NumMolecules>& num_molecules,
                           MolState<NumMolecules>& num_molecules_lower_bound,
                           MolState<NumMolecules>& num_molecules_upper_bound) const;
  /**
   * \brief Check all reaction rates bounds and update those as necessary.
   *
   * \param mol_state molecular state
   */
  void updateReactionRatesBounds(const MolState<NumMolecules> &mol_state,
                                 const osh::Real state_time);

  /**
   * \brief Check and update the reaction rates as necessary after a reaction
   * that happened in elementId.
   *
   * \param a_lower_bound lower bound
   * \param a_upper_bound upper_bound
   * \param mol_state molecular state
   * \param mol_state_element_updates updates.
   */
  void checkAndUpdateReactionRatesBounds(
      propensities_groups_t<kproc::PropensitiesPolicy::direct_event |
                            kproc::PropensitiesPolicy::without_next_event>
          &a_lower_bound,
      propensities_groups_t<kproc::PropensitiesPolicy::direct_event |
                            kproc::PropensitiesPolicy::with_next_event>
          &a_upper_bound,
      const MolState<NumMolecules> &mol_state, const Event &event,
      const std::vector<MolStateElementID> &mol_state_element_updates);

  /**
   * \brief Bounds generator for molecules.
   *
   * Bounds are defined based on the static member delta_rel_ as (nc * (1 -
   * delta_rel_), nc * (1 + delta_rel_)) where nc is the number of molecules. A
   * special treatment is applied for a small number of molecules. \param nc
   * number of molecules \param pPoolLB Lower bound on the number of molecules
   * (updated by the procedure). \param pPoolUB Upper bound on the number of
   * molecules (updated by the procedure).
   */
  static inline void applyBounds(NumMolecules nc,
                                 MolState<NumMolecules>& pPoolLB,
                                 MolState<NumMolecules>& pPoolUB,
                                 const MolStateElementID& elemID);

  MolState<NumMolecules> &pMolState;
  kproc::KProcState &pKProcState;

  // bounds on molecule populations
  MolState<NumMolecules> mol_state_lower_bound_;
  MolState<NumMolecules> mol_state_upper_bound_;

  // dependent reactions
  kproc::KProcState::DependenciesMap dependent_reactions_;

  // reaction propensity rates
  kproc::Propensities<NumMolecules, kproc::PropensitiesPolicy::direct_without_next_event>
      a_lower_bound_;
  kproc::Propensities<NumMolecules, kproc::PropensitiesPolicy::direct_with_next_event>
      a_upper_bound_;

  // uniform distribution
  std::uniform_real_distribution<double> uniform_;

  RNG &rng_;

  osh::Reals potential_on_vertices_;

  // counter
  osh::I64 extent{};

  // approximate relative distance of the molecules lower/upper bounds from the
  // actual molecules count.
  static constexpr osh::Real delta_rel_ = 0.05;
};

//--------------------------------------------------

// explicit instantiation declarations
extern template class RSSAOperator<std::mt19937, osh::I32>;
extern template class RSSAOperator<std::mt19937, osh::I64>;

extern template class RSSAOperator<steps::rng::RNG, osh::I32>;
extern template class RSSAOperator<steps::rng::RNG, osh::I64>;

//--------------------------------------------------

} // namespace dist
} // namespace steps
