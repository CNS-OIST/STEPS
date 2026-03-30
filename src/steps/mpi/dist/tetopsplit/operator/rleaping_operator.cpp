#include "rleaping_operator.hpp"

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_id.hpp"
#include "util/collections.hpp"
#include "util/profile/profiler_interface.hpp"
#include <cstddef>
#include <limits>
#include <random>
#include <stdexcept>
#include <string>

#undef MPI_Allreduce

namespace steps::dist {

RLeapingOperator::RLeapingOperator(MolState& mol_state,
                                   kproc::KProcState& kproc_state,
                                   rng::RNG& t_rng)
    : pMolState(mol_state)
    , pKProcState(kproc_state)
    , rng_(t_rng) {
    ArgErrLogIf(pKProcState.hasComplexReactions(),
                "The RSSA operator cannot be used if the model contains complexes with "
                "statesAsSpecies=False.");
    pKProcState.initPropensities(pPropensities);
    for (auto& group: pPropensities.groups()) {
        group.init(pMolState, pKProcState);
    }
    L_values.resize(pPropensities.groups().size(), 0);
    nb_L_steps.resize(pPropensities.groups().size(), 0);
}

std::map<std::string, double> RLeapingOperator::getDebugInfo() const {
    std::map<std::string, double> info;
    info["total_ssa_steps"] = total_ssa_steps;
    info["total_rleaping_steps"] = total_rleaping_steps;
    info["total_ssa_time"] = total_ssa_time;
    info["total_rleaping_time"] = total_rleaping_time;
    info["total_negative_population"] = total_negative_population;
    info["total_rleaping_overshoot_steps"] = total_rleaping_overshoot_steps;
    return info;
}

void RLeapingOperator::run(const osh::Real period, const osh::Real state_time) {
    Instrumentor::phase p("RLeapingOperator::run()");
    for (uint i = 0; i < pPropensities.groups().size(); ++i) {
        auto& group = pPropensities.groups()[i];
        auto& L = L_values[i];
        auto& L_steps = nb_L_steps[i];

        auto& outdated_kprocs = pKProcState.get_outdated_kprocs(i);
        group.update(pMolState, rng_, outdated_kprocs);
        outdated_kprocs.clear();

        osh::Real cumulative_dt{};
        bool useBinomial = true;
        osh::Real tau;
        osh::Real kahan = 0;
        uint nb_ssa_steps{0};
        bool toEnd = false;
        bool redo = false;
        osh::Real a0{0.0};
        while (true) {
            a0 = group.a0();
            if (useBinomial) {
                if (redo) {
                    redo = false;
                    toEnd = false;
                } else {
                    if (L_steps++ % L_compute_period == 0) {
                        L = group.computeL(pMolState, tolerance, theta);
                    }
                }
                if (L < ssa_threshold) {
                    // If L is too low run a number of SSA steps in a row, avoiding to recompute L
                    // too frequently
                    useBinomial = false;
                    nb_ssa_steps = 0;
                    L = 1;
                    if (group.a0() == 0) {
                        break;
                    }
                    tau = rng_.getExp(group.a0());
                    ++total_ssa_steps;
                    total_ssa_time += std::min(tau, period - cumulative_dt);
                } else if (group.a0() > 0) {
                    // Sample the time that corresponds to L reactions
                    std::gamma_distribution<> dist(L, 1.0 / group.a0());
                    tau = dist(rng_);
                    ++total_rleaping_steps;
                    if (cumulative_dt + tau > period) {
                        // If the time overshoots the end time, reduce L
                        ++total_rleaping_overshoot_steps;
                        std::poisson_distribution<> dist2((period - cumulative_dt) * group.a0());
                        L = dist2(rng_);
                        if (L == 0) {
                            break;
                        }
                        // Set tau to zero to avoid exiting early, set toEnd to true so that, after
                        // processing the L reactions, we exit
                        tau = 0;
                        toEnd = true;
                        total_rleaping_time += period - cumulative_dt;
                    } else {
                        total_rleaping_time += tau;
                    }
                } else {
                    // If no reactions to fire, just exit the loop
                    break;
                }
            } else {
                // Standard SSA step
                L = 1;
                if (group.a0() == 0) {
                    break;
                }
                tau = rng_.getExp(group.a0());
                ++total_ssa_steps;
                total_ssa_time += std::min(tau, period - cumulative_dt);
                if (++nb_ssa_steps > ssa_steps) {
                    useBinomial = true;
                }
            }
            if (cumulative_dt + tau > period) {
                break;
            }

            std::vector<MolStateElementID> updates;
            std::vector<osh::LO> dependencies;
            updates.reserve(L);
            dependencies.reserve(L);
            auto events = group.drawEvents(rng_, L);

            for (auto& [kid, nb]: events) {
                auto& upd = pKProcState.updateMolStateAndOccupancy(
                    pMolState, rng_, state_time + cumulative_dt, kid, nb);
                pKProcState.updateChargeFlow(kid, nb);
                const auto& deps = pKProcState.dependenciesFromEvent(kid);
                updates.insert(updates.end(), upd.begin(), upd.end());
                dependencies.insert(dependencies.end(), deps.begin(), deps.end());
            }

            bool has_negative = false;
            // Only check for negative population if we fired more than one reaction
            if (L > 0) {
                // Only species that were updated need to be checked
                for (auto mseid: updates) {
                    if (pMolState(mseid) < 0) {
                        ++total_negative_population;
                        has_negative = true;
                        break;
                    }
                }
            }
            if (has_negative) {
                // Apply the events reversed to restore the original state
                for (auto& [kid, nb]: events) {
                    pKProcState.updateMolStateAndOccupancy(
                        pMolState, rng_, state_time + cumulative_dt, kid, -nb);
                    pKProcState.updateChargeFlow(kid, -nb);
                }
                redo = true;
                L = L / 2;
                if (L < 1) {
                    L = 1;
                }
            } else {
                group.update(pMolState, rng_, dependencies);

                extent += L;
                if (toEnd) {
                    cumulative_dt = period;
                    break;
                } else {
                    // Use Crank-Nicolson type scheme to draw tau
                    if (L >= ssa_threshold) {
                        std::gamma_distribution<> dist(L, 2.0 / (a0 + group.a0()));
                        tau = dist(rng_);
                    }

                    // Use Kahan summation to avoid getting stuck when cumulative_dt >> tau
                    osh::Real tau2 = tau - kahan;
                    osh::Real sum = cumulative_dt + tau2;
                    kahan = (sum - cumulative_dt) - tau2;
                    cumulative_dt = sum;
                }
            }
        }
    }
    need_reset = false;
}

void RLeapingOperator::reset() {
    need_reset = true;
}

void RLeapingOperator::resetAndUpdateAll(const osh::Real state_time, const osh::Real /*max_time*/) {
    for (auto& group: pPropensities.groups()) {
        if (need_reset) {
            group.reset(pMolState, rng_, state_time);
        }
        group.update_all(pMolState, rng_, state_time);
    }
    std::fill(nb_L_steps.begin(), nb_L_steps.end(), 0);
    need_reset = false;
}

}  // namespace steps::dist
