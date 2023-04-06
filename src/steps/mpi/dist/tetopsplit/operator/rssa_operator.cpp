#include "rssa_operator.hpp"

#include <functional>
#include <vector>

#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/kproc/diffusions.hpp"
#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/reactions.hpp"

#undef MPI_Allreduce

namespace steps {
namespace dist {

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
RSSAOperator<RNG, NumMolecules>::RSSAOperator(MolState<NumMolecules>& mol_state,
                                              kproc::KProcState<NumMolecules>& k_proc_state,
                                              RNG& t_rng,
                                              osh::Reals potential_on_vertices)
    : pMolState(mol_state)
    , pKProcState(k_proc_state)
    , mol_state_lower_bound_(pMolState.species_per_elements(),
                             false,
                             pMolState.species_per_boundaries())
    , mol_state_upper_bound_(pMolState.species_per_elements(),
                             false,
                             pMolState.species_per_boundaries())
    , uniform_(0.0, 1.0)
    , rng_(t_rng)
    , potential_on_vertices_(potential_on_vertices) {
    static_assert(delta_rel_ >= 0.0 && delta_rel_ <= 0.5, "delta_rel_ out of bounds");
    pKProcState.initPropensities(a_lower_bound_);
    pKProcState.initPropensities(a_upper_bound_);
    k_proc_state.collateAllDependencies(dependent_reactions_);
}

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
template <typename Entity>
void RSSAOperator<RNG, NumMolecules>::updateNumMolsBounds(
    const EntityMolecules<Entity, NumMolecules>& num_molecules,
    MolState<NumMolecules>& num_molecules_lower_bound,
    MolState<NumMolecules>& num_molecules_upper_bound) const {
    for (auto entity: num_molecules.entities()) {
        for (auto species: num_molecules.species(entity)) {
            applyBounds(num_molecules(entity, species),
                        num_molecules_lower_bound,
                        num_molecules_upper_bound,
                        {entity, species});
            assert(num_molecules_lower_bound(entity, species) <=
                   num_molecules_upper_bound(entity, species));
        }
    }
}

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
void RSSAOperator<RNG, NumMolecules>::updateReactionRatesBounds(
    const MolState<NumMolecules>& mol_state, const osh::Real state_time) {
    updateNumMolsBounds(mol_state.moleculesOnElements(),
                        mol_state_lower_bound_,
                        mol_state_upper_bound_);
    updateNumMolsBounds(mol_state.moleculesOnPatchBoundaries(),
                        mol_state_lower_bound_,
                        mol_state_upper_bound_);

    for (size_t k = 0; k < a_upper_bound_.groups().size(); k++) {
        a_upper_bound_.groups()[k].update_all(mol_state_upper_bound_, rng_, state_time);
        a_lower_bound_.groups()[k].update_all(mol_state_lower_bound_, rng_, state_time);
    }
}

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
inline void RSSAOperator<RNG, NumMolecules>::applyBounds(const NumMolecules nc,
                                                         MolState<NumMolecules>& low,
                                                         MolState<NumMolecules>& up,
                                                         const MolStateElementID& elemID) {
    const auto nc_real = static_cast<osh::Real>(nc);

    if (nc_real > 3.0 / delta_rel_) {
        low.assign(elemID, static_cast<NumMolecules>(nc_real * (1 - delta_rel_)));
        up.assign(elemID, static_cast<NumMolecules>(nc_real * (1 + delta_rel_)));
    } else if (nc > 6) {
        // nc >= 7 && nc <= 3.0 / delta_rel_
        low.assign(elemID, nc - 3);
        up.assign(elemID, nc + 3);
    } else {
        low.assign(elemID, nc);
        up.assign(elemID, nc);
    }
}

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
void RSSAOperator<RNG, NumMolecules>::checkAndUpdateReactionRatesBounds(
    propensities_groups_t<kproc::PropensitiesPolicy::direct_without_next_event>& a_lower_bound,
    propensities_groups_t<kproc::PropensitiesPolicy::direct_with_next_event>& a_upper_bound,
    const MolState<NumMolecules>& mol_state,
    const Event& event,
    const std::vector<MolStateElementID>& mol_state_element_updates) {
    //  indices of reactions to recompute
    std::set<kproc::KProcID> reactions_to_recompute;
    for (auto el: mol_state_element_updates) {
        NumMolecules m = mol_state(el);
        const NumMolecules& upper = mol_state_upper_bound_(el);
        const NumMolecules& lower = mol_state_lower_bound_(el);
        if (m > upper || m < lower) {
            // recenter
            applyBounds(m, mol_state_lower_bound_, mol_state_upper_bound_, el);
            assert(lower <= upper);
            std::vector<kproc::KProcID>& dep_reacs = dependent_reactions_[el];
            reactions_to_recompute.insert(dep_reacs.begin(), dep_reacs.end());
        }
    }

    a_lower_bound.template update<std::set<kproc::KProcID>>(mol_state_lower_bound_,
                                                            rng_,
                                                            event,
                                                            reactions_to_recompute);
    a_upper_bound.template update<std::set<kproc::KProcID>>(mol_state_upper_bound_,
                                                            rng_,
                                                            event,
                                                            reactions_to_recompute);
}

//---------------------------------------------------------

template <typename RNG, typename NumMolecules>
osh::Real RSSAOperator<RNG, NumMolecules>::run(const osh::Real period, const osh::Real state_time) {
    updateReactionRatesBounds(pMolState, state_time);  // setup the propensity lower/upper bounds.
    osh::Real slack{-period};
    for (size_t k = 0; k < a_upper_bound_.groups().size(); ++k) {
        osh::Real cumulative_dt(0.0);
        while (true) {
            bool isRejected(true);
            unsigned kproc_id_data;
            do {
                std::pair<osh::Real, kproc::KProcID> next_event =
                    a_upper_bound_.groups()[k].drawEvent(rng_, cumulative_dt);
                cumulative_dt = next_event.first;
                kproc::KProcID& reaction_candidate_id = next_event.second;
                kproc_id_data = reaction_candidate_id.data();
                if (std::isinf(cumulative_dt)) {
                    break;
                }
                auto random_selector = uniform_(rng_) * a_upper_bound_[reaction_candidate_id];
                if (random_selector <= a_lower_bound_[reaction_candidate_id]) {
                    isRejected = false;
                } else {
                    // compute the 'true' propensity
                    osh::Real rate = pKProcState.propensityFun()(reaction_candidate_id, pMolState);
                    if (rate < std::numeric_limits<osh::Real>::epsilon()) {
                        throw std::logic_error(
                            "RSSA: propensity rate of the candidate reaction is zero.");
                    }
                    if (random_selector <= rate) {
                        isRejected = false;
                    }
                }

            } while (isRejected);
            kproc::KProcID reaction_id{kproc_id_data};
            if (cumulative_dt > period) {
                break;
            }
            slack = std::max(slack, cumulative_dt - period);
            if (reaction_id.type() == kproc::KProcType::GHKSReac) {
                pKProcState.updateGHKChargeFlow(reaction_id.id());
            }
            auto elementsUpdated = pKProcState.updateMolStateAndOccupancy(
                pMolState, state_time + cumulative_dt, reaction_id);
            checkAndUpdateReactionRatesBounds(a_lower_bound_.groups()[k],
                                              a_upper_bound_.groups()[k],
                                              pMolState,
                                              {cumulative_dt, reaction_id},
                                              elementsUpdated);
            extent += 1;
        }
    }
    return slack;
}

//---------------------------------------------------------

// explicit instantiation definitions
template class RSSAOperator<std::mt19937, osh::I32>;
template class RSSAOperator<std::mt19937, osh::I64>;
template class RSSAOperator<steps::rng::RNG, osh::I32>;
template class RSSAOperator<steps::rng::RNG, osh::I64>;
//---------------------------------------------------------

} // namespace dist
} // namespace steps
