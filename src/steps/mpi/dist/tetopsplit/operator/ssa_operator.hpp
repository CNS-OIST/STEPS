#pragma once

#include <random>

#include "../kproc/fwd.hpp"
#include "geom/dist/fwd.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/propensities.hpp"
#include "rng/rng.hpp"

namespace steps::dist {

/**
 * \brief SSA operator.
 *
 * This class implements the SSA operator in the OpSplit solution.
 *
 * The core function of this class is run(), which partially corresponds to the
 * SSA component the simulation core loop tetexact::Tetexact::run() /
 * tetopsplit::TetOpSplit::run() in STEPS.
 */
template <NextEventSearchMethod SearchMethod>
class SSAOperator {
  public:
    SSAOperator(MolState& mol_state,
                kproc::KProcState& kproc_state,
                rng::RNG& t_rng,
                osh::Reals potential_on_vertices);

    /// \return the number of kinetic events that arose on the current process
    inline osh::I64 getExtent() const noexcept {
        return extent;
    }

    /**
     * \brief Execute the operator
     *
     * \param period How long the operator should be executed
     * \param state_time The start state time when the operator is executed
     * \return osh::Real The slack between the cumulated dts of all ssa
     * events and the period
     */
    osh::Real run(const osh::Real period, const osh::Real state_time);

    /**
     * \brief Select the next kinetic process(reaction, surface reaction) in the
     * SSA.
     * \return The id of the next kinetic process
     */
    template <class Generator>
    kproc::KProcID getNextEvent(Generator& rng) const;

    /**
     * \brief reset the required propensity data stored in the operator.
     *
     * This function set the need_reset flag in the operator to true.
     * At the start of run() the operator update itself if the need_reset
     * flag is true, then set it to false.
     *
     * This delay reset is needed as the call of simulation reset is before
     * molecule/propensity changes, but the reset of this operator should be
     * performed after the molecule/propensity changes.
     */
    void reset();

    /**
     * \brief 1. Reset the group data structure (group of propensities)
     *        2. Set the maximum time to be considered as valid in the SSA (used in GB optimization)
     *        3. Update ALL propensities
     */
    void resetAndUpdateAll(const osh::Real state_time, const osh::Real max_time);

  private:
    MolState& pMolState;
    kproc::KProcState& pKProcState;
    rng::RNG& rng_;
    osh::Reals potential_on_vertices_;

    osh::I64 extent{};
    kproc::Propensities<kproc::PropensitiesPolicy::get<SearchMethod>() |
                        kproc::PropensitiesPolicy::with_next_event>
        pPropensities;

    bool need_reset{true};
};

// explicit template instantiation declarations
extern template class SSAOperator<NextEventSearchMethod::Direct>;
extern template class SSAOperator<NextEventSearchMethod::GibsonBruck>;

}  // namespace steps::dist
