#include "ssa_operator.hpp"

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_state.hpp"
#include "util/profile/profiler_interface.hpp"
#include <type_traits>

#undef MPI_Allreduce

namespace steps::dist {

template <NextEventSearchMethod SearchMethod>
SSAOperator<SearchMethod>::SSAOperator(MolState& mol_state,
                                       kproc::KProcState& kproc_state,
                                       rng::RNG& t_rng)
    : pMolState(mol_state)
    , pKProcState(kproc_state)
    , rng_(t_rng) {
    pKProcState.initPropensities(pPropensities);
}

template <NextEventSearchMethod SearchMethod>
osh::Real SSAOperator<SearchMethod>::run(const osh::Real period, const osh::Real state_time) {
    Instrumentor::phase p("SSAOperator::run()");
    osh::Real slack{-period};
    for (uint i = 0; i < pPropensities.groups().size(); ++i) {
        auto& group = pPropensities.groups()[i];

        auto& outdated_kprocs = pKProcState.get_outdated_kprocs(i);
        group.update_outdated(outdated_kprocs, pMolState, rng_, state_time);
        outdated_kprocs.clear();

        osh::Real cumulative_dt{};
        while (true) {
            const kproc::Event& event = group.drawEvent(rng_, state_time + cumulative_dt);
            if (event.first > (state_time + period)) {
                break;
            }
            cumulative_dt = (event.first - state_time);
            slack = std::max(slack, cumulative_dt - period);
            pKProcState.updateMolStateAndOccupancy(pMolState, rng_, event.first, event.second);
            pKProcState.updateChargeFlow(event.second);
            const kproc::KProcDeps& dependencies = pKProcState.dependenciesFromEvent(event.second);
            group.template update<kproc::KProcDeps>(pMolState, rng_, event, dependencies);
            extent += 1;
        }
    }
    need_reset = false;
    return slack;
}

template <NextEventSearchMethod SearchMethod>
void SSAOperator<SearchMethod>::reset() {
    need_reset = true;
}

template <NextEventSearchMethod SearchMethod>
void SSAOperator<SearchMethod>::resetAndUpdateAll(const osh::Real state_time,
                                                  const osh::Real max_time) {
    for (auto& group: pPropensities.groups()) {
        if (need_reset) {
            group.reset(pMolState, rng_, state_time);
        }
        if constexpr (SearchMethod == NextEventSearchMethod::GibsonBruck) {
            group.updateMaxTime(max_time);
        }
        group.update_all(pMolState, rng_, state_time);
    }
    need_reset = false;
}

// explicit template instantiation definitions
template class SSAOperator<NextEventSearchMethod::Direct>;
template class SSAOperator<NextEventSearchMethod::GibsonBruck>;

}  // namespace steps::dist
