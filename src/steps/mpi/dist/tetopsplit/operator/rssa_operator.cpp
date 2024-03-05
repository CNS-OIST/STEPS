#include "rssa_operator.hpp"

#include <functional>
#include <vector>

#include "mpi/dist/tetopsplit/definition/reacdef.hpp"
#include "mpi/dist/tetopsplit/kproc/diffusions.hpp"
#include "mpi/dist/tetopsplit/kproc/fwd.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/reactions.hpp"

#undef MPI_Allreduce

namespace steps::dist {

//---------------------------------------------------------

RSSAOperator::RSSAOperator(MolState& mol_state,
                           kproc::KProcState& k_proc_state,
                           rng::RNG& t_rng,
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
}

//---------------------------------------------------------

template <typename Entity>
void RSSAOperator::updateNumMolsBounds(const EntityMolecules<Entity>& num_molecules,
                                       MolState& num_molecules_lower_bound,
                                       MolState& num_molecules_upper_bound) const {
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

void RSSAOperator::updateReactionRatesBounds(const MolState& mol_state,
                                             const osh::Real state_time) {
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

inline void RSSAOperator::applyBounds(const molecules_t nc,
                                      MolState& low,
                                      MolState& up,
                                      const MolStateElementID& elemID) {
    const auto nc_real = static_cast<osh::Real>(nc);

    if (nc_real > 3.0 / delta_rel_) {
        low.assign(elemID, static_cast<molecules_t>(nc_real * (1 - delta_rel_)));
        up.assign(elemID, static_cast<molecules_t>(nc_real * (1 + delta_rel_)));
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

void RSSAOperator::checkAndUpdateReactionRatesBounds(
    propensities_groups_t<kproc::PropensitiesPolicy::direct_without_next_event>& a_lower_bound,
    propensities_groups_t<kproc::PropensitiesPolicy::direct_with_next_event>& a_upper_bound,
    const MolState& mol_state,
    const Event& event,
    const std::vector<MolStateElementID>& mol_state_element_updates) {
    //  indices of reactions to recompute
    std::set<kproc::KProcID> reactions_to_recompute;
    for (auto el: mol_state_element_updates) {
        molecules_t m = mol_state(el);
        const molecules_t& upper = mol_state_upper_bound_(el);
        const molecules_t& lower = mol_state_lower_bound_(el);
        if (m > upper || m < lower) {
            // recenter
            applyBounds(m, mol_state_lower_bound_, mol_state_upper_bound_, el);
            assert(lower <= upper);

            auto species = std::get<1>(el);
            std::visit(
                [&](auto& entity) {
                    using T = std::decay_t<decltype(entity)>;

                    if constexpr (std::is_same_v<T, mesh::tetrahedron_id_t>) {
                        const auto idx = pMolState.moleculesOnElements().ab(entity, species);
                        for (auto elem: pKProcState.get_dependency_map_elems()[idx]) {
                            reactions_to_recompute.emplace(elem);
                        }
                    } else if constexpr (std::is_same_v<T, mesh::triangle_id_t>) {
                        const auto idx = pMolState.moleculesOnPatchBoundaries().ab(entity, species);
                        for (auto bound: pKProcState.get_dependency_map_bnds()[idx]) {
                            reactions_to_recompute.emplace(bound);
                        }
                    } else {
                        static_assert(util::always_false_v<T>, "unmanaged entity type");
                    }
                },
                std::get<0>(el));
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

osh::Real RSSAOperator::run(const osh::Real period, const osh::Real state_time) {
    updateReactionRatesBounds(pMolState, state_time);  // setup the propensity lower/upper bounds.
    osh::Real slack{-period};
    for (size_t k = 0; k < a_upper_bound_.groups().size(); ++k) {
        osh::Real cumulative_dt(0.0);
        while (true) {
            bool isRejected(true);
            unsigned kproc_id_data;
            do {
                const std::pair<osh::Real, kproc::KProcID>& next_event =
                    a_upper_bound_.groups()[k].drawEvent(rng_, cumulative_dt);
                cumulative_dt = next_event.first;
                const kproc::KProcID& reaction_candidate_id = next_event.second;
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

}  // namespace steps::dist
