#pragma once

#include <random>

#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/propensities.hpp"
#include "mpi/dist/tetopsplit/kproc/reactions.hpp"
#include "rng/rng.hpp"
#include "util/collections.hpp"

namespace steps::dist {

/**
 * \brief RSSA operator.
 *
 * This class implements the RSSA operator in the OpSplit solution.
 *
 * The core function of this class is run(), which partially corresponds to the
 * SSA component the simulation core loop tetexact::Tetexact::run() /
 * tetopsplit::TetOpSplit::run() in STEPS.
 */
class RSSAOperator {
  public:
    template <int Policy>
    using propensities_groups_t = kproc::PropensitiesGroup<Policy>;
    using Event = std::pair<osh::Real, kproc::KProcID>;
    /**
     * Constructor
     * \param mol_state molecular state
     * \param dv molecules transferred to neighbours
     * \param kproc_state kinetic processes state
     */
    RSSAOperator(MolState& mol_state,
                 kproc::KProcState& kproc_state,
                 rng::RNG& t_rng,
                 osh::Reals potential_on_vertices);

    /// \return the number of kinetic events that arose on the current process
    inline osh::I64 getExtent() const noexcept {
        return extent;
    }

    osh::Real run(osh::Real period, osh::Real state_time);

    /* TODO : need further implementation and testing for the new GB optimizations */
    /**
     * \brief reset the required propensity data stored in the operator
     */
    void reset() {
        // make it do nothing now to pass the CI.
    }

    /* TODO : need further implementation and testing for the new GB optimizations */
    /**
     * \brief 1. Reset the group data structure (group of propensities)
     *        2. Set the maximum time to be considered as valid in the SSA (used in GB optimization)
     *        3. Update ALL propensities
     */
    void resetAndUpdateAll(const osh::Real /*state_time*/, const osh::Real /*max_time*/) {
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
    void updateNumMolsBounds(const EntityMolecules<Entity>& num_molecules,
                             MolState& num_molecules_lower_bound,
                             MolState& num_molecules_upper_bound) const;
    /**
     * \brief Check all reaction rates bounds and update those as necessary.
     *
     * \param mol_state molecular state
     */
    void updateReactionRatesBounds(const MolState& mol_state, const osh::Real state_time);

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
        propensities_groups_t<kproc::PropensitiesPolicy::direct_without_next_event>& a_lower_bound,
        propensities_groups_t<kproc::PropensitiesPolicy::direct_with_next_event>& a_upper_bound,
        const MolState& mol_state,
        const Event& event,
        const std::vector<MolStateElementID>& mol_state_element_updates);

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
    static inline void applyBounds(molecules_t nc,
                                   MolState& pPoolLB,
                                   MolState& pPoolUB,
                                   const MolStateElementID& elemID);

    MolState& pMolState;
    kproc::KProcState& pKProcState;

    // bounds on molecule populations
    MolState mol_state_lower_bound_;
    MolState mol_state_upper_bound_;

    // reaction propensity rates
    kproc::Propensities<kproc::PropensitiesPolicy::direct_without_next_event> a_lower_bound_;
    kproc::Propensities<kproc::PropensitiesPolicy::direct_with_next_event> a_upper_bound_;

    // uniform distribution
    std::uniform_real_distribution<double> uniform_;

    rng::RNG& rng_;

    osh::Reals potential_on_vertices_;

    // counter
    osh::I64 extent{};

    // approximate relative distance of the molecules lower/upper bounds from the
    // actual molecules count.
    static constexpr osh::Real delta_rel_ = 0.05;
};

}  // namespace steps::dist
