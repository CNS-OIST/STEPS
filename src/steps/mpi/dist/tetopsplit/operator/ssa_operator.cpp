#include "ssa_operator.hpp"

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_state.hpp"
#include "util/profile/profiler_interface.h"

#undef MPI_Allreduce

namespace steps {
namespace dist {

template <typename RNG, typename NumMolecules, NextEventSearchMethod SearchMethod>
SSAOperator<RNG, NumMolecules, SearchMethod>::SSAOperator(
    MolState<NumMolecules>& mol_state,
    kproc::KProcState<NumMolecules>& kproc_state,
    RNG& t_rng,
    osh::Reals potential_on_vertices)
    : pMolState(mol_state)
    , pKProcState(kproc_state)
    , rng_(t_rng)
    , potential_on_vertices_(potential_on_vertices) {
    pKProcState.initPropensities(pPropensities);
}

template <typename RNG, typename NumMolecules, NextEventSearchMethod SearchMethod>
osh::Real SSAOperator<RNG, NumMolecules, SearchMethod>::run(const osh::Real period,
                                                            const osh::Real state_time) {
    Instrumentor::phase p("SSAOperator::run()");
    osh::Real slack{-period};
    for (auto& group: pPropensities.groups()) {
        group.template update_outdated<RNG>(pMolState, rng_, state_time);
        osh::Real cumulative_dt{};
        while (true) {
            const kproc::Event& event = group.drawEvent(rng_, state_time + cumulative_dt);
            if (event.first > (state_time + period)) {
                break;
            }
            cumulative_dt = (event.first - state_time);
            slack = std::max(slack, cumulative_dt - period);
            pKProcState.updateMolStateAndOccupancy(pMolState,
                                                   event.first,
                                                   event.second);
            if (event.second.type() == kproc::KProcType::GHKSReac) {
                pKProcState.updateGHKChargeFlow(event.second.id());
            }
            const kproc::KProcDeps& dependencies = pKProcState.dependenciesFromEvent(event.second);
            group.template update<kproc::KProcDeps, RNG>(pMolState, rng_, event, dependencies);
            extent += 1;
        }
    }
    need_reset = false;
    return slack;
}

template <typename RNG, typename NumMolecules, NextEventSearchMethod SearchMethod>
void SSAOperator<RNG, NumMolecules, SearchMethod>::reset() {
    need_reset = true;
}

template <typename RNG, typename NumMolecules, NextEventSearchMethod SearchMethod>
void SSAOperator<RNG, NumMolecules, SearchMethod>::resetAndUpdateAll(const osh::Real state_time,
                                                                     const osh::Real max_time) {
    for (auto& group: pPropensities.groups()) {
        if (need_reset) {
            group.template reset<RNG>(pMolState, rng_, state_time);
        }
        group.updateMaxTime(max_time);
        group.template update_all<RNG>(pMolState, rng_, state_time);
    }
    need_reset = false;
}

// explicit template instantiation definitions
template class SSAOperator<std::mt19937, osh::I32,
                           NextEventSearchMethod::Direct>;
template class SSAOperator<std::mt19937, osh::I64,
                           NextEventSearchMethod::Direct>;
template class SSAOperator<std::mt19937, osh::I32,
                           NextEventSearchMethod::GibsonBruck>;
template class SSAOperator<std::mt19937, osh::I64,
                           NextEventSearchMethod::GibsonBruck>;

template class SSAOperator<steps::rng::RNG, osh::I32,
                           NextEventSearchMethod::Direct>;
template class SSAOperator<steps::rng::RNG, osh::I64,
                           NextEventSearchMethod::Direct>;
template class SSAOperator<steps::rng::RNG, osh::I32,
                           NextEventSearchMethod::GibsonBruck>;
template class SSAOperator<steps::rng::RNG, osh::I64,
                           NextEventSearchMethod::GibsonBruck>;

} // namespace dist
} // namespace steps
