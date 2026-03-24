#pragma once


#include "../kproc/fwd.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/propensities.hpp"
#include "rng/rng.hpp"

namespace steps::dist {

/**
 * \brief Reaction leaping operator
 *
 * This class implements the SSA R-leaping operator in the OpSplit solution.
 *
 * A. Auger, P. Chatelain, and P. Koumoutsakos. R-leaping: Accelerating the stochastic simulation
 * algorithm by reaction leaps. The Journal of chemical physics, 125(8), 2006.
 */
class RLeapingOperator {
  public:
    RLeapingOperator(MolState& mol_state, kproc::KProcState& kproc_state, rng::RNG& t_rng);

    /// \return the number of kinetic events that arose on the current process
    inline osh::I64 getExtent() const noexcept {
        return extent;
    }

    /**
     * \brief Get debugging information
     *
     * Return the following values:
     *   - total_ssa_steps: The total number of standard SSA steps (steps that involve only one
     *                      reaction)
     *   - total_rleaping_steps: The total number of RLeaping steps (steps that involve several
     *                           reactions)
     *   - total_ssa_time: The total biological time covered by standard SSA steps
     *   - total_rleaping_time: The total biological time covered by RLeaping steps
     *   - total_negative_population: The total number of RLeaping steps that resulted in negative
     *                                species populations
     *   - total_rleaping_overshoot_steps: The total number of RLeaping steps that were reduced
     * because they would have resulted in a jump further than the required end time
     *
     * \return a map from value name to value
     */
    std::map<std::string, double> getDebugInfo() const;

    /**
     * \brief Execute the operator
     *
     * \param period How long the operator should be executed
     * \param state_time The start state time when the operator is executed
     */
    void run(const osh::Real period, const osh::Real state_time);

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
     * \brief Reset the group data structure (group of propensities) and update all propensities
     */
    void resetAndUpdateAll(const osh::Real state_time, const osh::Real max_time);

    //////////////////////////////////////
    // Getters / setters for parameters //
    //////////////////////////////////////

    inline uint getSSAThreshold() const noexcept {
        return ssa_threshold;
    }
    inline void setSSAThreshold(uint thresh) noexcept {
        ssa_threshold = thresh;
    }

    inline osh::Real getTolerance() const noexcept {
        return tolerance;
    }
    inline void setTolerance(osh::Real _tolerance) {
        ArgErrLogIf(_tolerance <= 0, "Tolerance should be strictly positive");
        tolerance = _tolerance;
    }

    inline osh::Real getTheta() const noexcept {
        return theta;
    }
    inline void setTheta(osh::Real _theta) {
        ArgErrLogIf(_theta < 0, "Theta should be positive");
        theta = _theta;
    }

    inline uint getSSASteps() const noexcept {
        return ssa_steps;
    }
    inline void setSSASteps(uint steps) noexcept {
        ssa_steps = steps;
    }

    inline uint getLComputePeriod() const noexcept {
        return L_compute_period;
    }
    inline void setLComputePeriod(uint period) noexcept {
        L_compute_period = period;
    }

  private:
    MolState& pMolState;
    kproc::KProcState& pKProcState;
    rng::RNG& rng_;

    osh::I64 extent{};

    // Total number of SSA steps (single steps)
    osh::I64 total_ssa_steps{0};
    // Total number of RLeaping steps (multiple SSA steps at once)
    osh::I64 total_rleaping_steps{0};
    // Total simulation time covered by SSA steps
    osh::Real total_ssa_time{0};
    // Total simulation time covered by RLeaping steps
    osh::Real total_rleaping_time{0};
    // Total number of times that an RLeaping step resulted in negative populations
    osh::I64 total_negative_population{0};
    // Total number of RLeaping steps that lead to overshooting the end period
    osh::I64 total_rleaping_overshoot_steps{0};

    kproc::Propensities<kproc::PropensitiesPolicy::get<NextEventSearchMethod::RLeaping>() |
                        kproc::PropensitiesPolicy::rleaping_without_next_event>
        pPropensities;

    // L values for each propensity group
    std::vector<long long> L_values;
    // Number of times L was used without recomputation
    std::vector<uint> nb_L_steps;

    // Tolerance parameter for computing L (the number of simultaneous reactions)
    double tolerance{0.05};
    // Parameter for controlling the apparition of negative species
    double theta{0.1};
    // Number of SSA steps to run in a row when L is too low
    uint ssa_steps{100};
    // The minimum number of simultaneous reactions for a reaction leaping step
    uint ssa_threshold{50};
    // The period at which L is computed. 1 means L is always recomputed.
    uint L_compute_period{10};

    bool need_reset{true};
};

}  // namespace steps::dist
