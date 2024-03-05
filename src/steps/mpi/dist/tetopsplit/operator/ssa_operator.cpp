#include "ssa_operator.hpp"

#include "../kproc/diffusions.hpp"
#include "../kproc/kproc_state.hpp"
#include "util/profile/profiler_interface.hpp"

#undef MPI_Allreduce

namespace steps::dist {

template <NextEventSearchMethod SearchMethod>
SSAOperator<SearchMethod>::SSAOperator(MolState& mol_state,
                                       kproc::KProcState& kproc_state,
                                       rng::RNG& t_rng,
                                       osh::Reals potential_on_vertices)
    : pMolState(mol_state)
    , pKProcState(kproc_state)
    , rng_(t_rng)
    , potential_on_vertices_(potential_on_vertices) {
    pKProcState.initPropensities(pPropensities);
}

template <NextEventSearchMethod SearchMethod>
osh::Real SSAOperator<SearchMethod>::run(const osh::Real period, const osh::Real state_time) {
    Instrumentor::phase p("SSAOperator::run()");
    osh::Real slack{-period};
    for (auto& group: pPropensities.groups()) {
        group.update_outdated(pMolState, rng_, state_time);
        osh::Real cumulative_dt{};
        while (true) {
            const kproc::Event& event = group.drawEvent(rng_, state_time + cumulative_dt);
            if (event.first > (state_time + period)) {
                break;
            }
            cumulative_dt = (event.first - state_time);
            slack = std::max(slack, cumulative_dt - period);
            pKProcState.updateMolStateAndOccupancy(pMolState, event.first, event.second);
            if (event.second.type() == kproc::KProcType::GHKSReac) {
                pKProcState.updateGHKChargeFlow(event.second.id());
            }
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
        group.updateMaxTime(max_time);
        group.update_all(pMolState, rng_, state_time);
    }
    need_reset = false;
}

// explicit template instantiation definitions
template class SSAOperator<NextEventSearchMethod::Direct>;
template class SSAOperator<NextEventSearchMethod::GibsonBruck>;

}  // namespace steps::dist
