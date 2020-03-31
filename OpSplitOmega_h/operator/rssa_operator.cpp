#include "rssa_operator.hpp"

#include <functional>
#include <vector>

#include "../kproc/diffusions.hpp"
#include "../kproc/reactions.hpp"
#include "mpitools.hpp"
#include "opsplit/reacdef.hpp"

#undef MPI_Allreduce

namespace zee {


//---------------------------------------------------------
template <typename RNG>
template <osh::Int Dim>
RSSAOperator<RNG>::RSSAOperator(MolState& mol_state,
                                DiffusionVariables& dv,
                                KProcState& kproc_state,
                                const OmegaHMesh<Dim>& mesh,
                                RNG& t_rng)
    : pMolState(mol_state)
    , pDiffusionVariables(dv)
    , pKProcState(kproc_state)
    , owned_elements_(mesh.getOwnedElems())
    , mol_state_lower_bound_(pMolState.species_per_elements())
    , mol_state_upper_bound_(pMolState.species_per_elements())
    , a_lower_bound_(kproc_state.reactions().size(), 0.0)
    , uniform_(0.0, 1.0)
    , rng_(t_rng) {
    if (kproc_state.surfaceReactions().size() > 0) {
        throw std::invalid_argument("RSSA operator cannot be used with surface reactions.");
    }
    static_assert(delta_rel_ >= 0.0 && delta_rel_ <= 0.5, "delta_rel_ out of bounds");

    // set up the dependent reactions map
    for (size_t k = 0; k < kproc_state.reactions().size(); ++k) {
        const auto& reacDef = kproc_state.reactions().getReacDef(k);
        const auto& lhs = reacDef.getPoolChangeArray(Reacdef::PoolChangeArrayType::LHS);
        const auto ownedPoint = kproc_state.reactions().getOwnerPoint(k);
        for (container::specie_id specie_id{}; specie_id.get() < static_cast<osh::LO>(lhs.size());
             ++specie_id) {
            if (lhs[static_cast<size_t>(specie_id.get())] != 0) {
                dependent_reactions_[std::make_pair(ownedPoint, specie_id)].insert(k);
            }
        }
    }

    // initialize bounds on firing rates, based on lower/upper bounds on mols
    for (auto e: owned_elements_) {
        for (container::specie_id s(0); s.get() < pMolState.numSpecies(e); ++s) {
            applyBounds(mol_state(e, s),
                        mol_state_lower_bound_(e, s),
                        mol_state_upper_bound_(e, s));
            assert(mol_state_lower_bound_(e, s) <= mol_state_upper_bound_(e, s));
        }
    }
    for (size_t index{}; index < kproc_state.reactions().size(); ++index) {
        a_lower_bound_[index] = kproc_state.reactions().computeRate(mol_state_lower_bound_, index);
    }
    std::vector<PetscScalar> a_upper_bound_temp(kproc_state.reactions().size());
    for (size_t index{}; index < kproc_state.reactions().size(); ++index) {
        a_upper_bound_temp[index] = kproc_state.reactions().computeRate(mol_state_upper_bound_,
                                                                        index);
    }
    a_upper_bound_.init(a_upper_bound_temp);
}

//---------------------------------------------------------

template <typename RNG>
inline void RSSAOperator<RNG>::applyBounds(const osh::LO nc, osh::LO& low, osh::LO& up) {
    if (nc > 3.0 / delta_rel_) {
        low = static_cast<osh::LO>(nc * (1 - delta_rel_));
        up = static_cast<osh::LO>(nc * (1 + delta_rel_));
    } else if (nc > 3) {
        // nc >= 4 && nc <= 3.0 / delta_rel_
        low = nc - 3;
        up = nc + 3;
    } else if (nc > 0) {
        // nc = 1 or 2 or 3
        low = 1;
        up = 2 * nc;
    } else {
        // nc = 0
        low = 0;
        up = 0;
    }
}

//---------------------------------------------------------

template <typename RNG>
void RSSAOperator<RNG>::checkAndUpdateReactionRatesBounds(
    const MolState& mol_state,
    const mesh::element_id elementId,
    const boost::optional<std::vector<PetscInt>>& reaction_indices) {
    //  indices of reactions to recompute
    std::set<size_t> reactions_to_recompute;
    for (container::specie_id specie{}; specie.get() < mol_state.numSpecies(elementId); ++specie) {
        if (!reaction_indices || (*reaction_indices)[static_cast<size_t>(specie.get())] != 0) {
            osh::LO m = mol_state(elementId, specie);
            osh::LO& upper = mol_state_upper_bound_(elementId, specie);
            osh::LO& lower = mol_state_lower_bound_(elementId, specie);
            if (m > upper || m < lower) {
                // recenter
                applyBounds(m, lower, upper);
                assert(lower <= upper);
                std::set<size_t>& dep_reacs =
                    dependent_reactions_[std::make_pair(elementId, specie)];
                reactions_to_recompute.insert(dep_reacs.begin(), dep_reacs.end());
            }
        }
    }
    for (auto& idx: reactions_to_recompute) {
        a_lower_bound_[idx] = pKProcState.reactions().computeRate(mol_state_lower_bound_, idx);
        a_upper_bound_.update(idx,
                              pKProcState.reactions().computeRate(mol_state_upper_bound_, idx));
    }
}

//---------------------------------------------------------

template <typename RNG>
void RSSAOperator<RNG>::checkAndUpdateReactionRatesBounds(const MolState& mol_state) {
    // check whether we need to update the mol state bounds
    for (const auto element: owned_elements_) {
        checkAndUpdateReactionRatesBounds(mol_state, element);
    }
}

//---------------------------------------------------------

// TODO(smelchio): alternative to avoid copy pate
template <typename RNG>
PetscInt64 RSSAOperator<RNG>::getExtent() const {
    PetscInt64 global_ext;
    MPI_Allreduce(&extent, &global_ext, 1, MPI_INT64_T, MPI_SUM, MPI_COMM_WORLD);
    return global_ext;
}

//---------------------------------------------------------

template <typename RNG>
void RSSAOperator<RNG>::run(PetscScalar period) {
    checkAndUpdateReactionRatesBounds(pMolState);  // setup the propensity lower/upper bounds.
    PetscScalar cumulative_dt(0.0);

    while (true) {
        if (a_upper_bound_.sum() < -std::numeric_limits<PetscScalar>::epsilon()) {
            throw std::invalid_argument("Negative propensity encountered.");
        }
        if (a_upper_bound_.sum() < std::numeric_limits<PetscScalar>::epsilon()) {
            // reaction can't happen, exit
            break;
        }

        bool isRejected(true);
        double erlangFactor(1.0);
        size_t reaction_candidate_index;

        do {
            reaction_candidate_index = singleOutReaction();

            auto random_selector = uniform_(rng_) * a_upper_bound_[reaction_candidate_index];
            if (random_selector <= a_lower_bound_[reaction_candidate_index]) {
                isRejected = false;
            } else {
                // compute the 'true' propensity
                const PetscScalar rate =
                    pKProcState.reactions().computeRate(pMolState, reaction_candidate_index);

                if (rate < std::numeric_limits<PetscScalar>::epsilon()) {
                    throw std::logic_error(
                        "RSSA: propensity rate of the candidate reaction is zero.");
                }
                if (random_selector <= rate) {
                    isRejected = false;
                }
            }
            erlangFactor *= uniform_(rng_);
        } while (isRejected);

        // const reaction_ptr& next = reacs[reaction_candidate_index];
        double dt = -1.0 / a_upper_bound_.sum() * std::log(erlangFactor);
        if (cumulative_dt + dt > period) {
            break;
        }
        cumulative_dt += dt;
        const auto ownerPoint = pKProcState.reactions().getOwnerPoint(reaction_candidate_index);
        // SMR: pDiffusionVariables might output a struct to have only one access =>
        // gathering all data in 1 object encapsulating two Kokkos arrays
        pKProcState.updateOccupancy(pDiffusionVariables,
                                    pMolState,
                                    cumulative_dt,
                                    std::make_pair(KProcType::Reac, reaction_candidate_index));
        pKProcState.reactions().apply(pMolState, reaction_candidate_index);
        const auto& reacDef = pKProcState.reactions().getReacDef(reaction_candidate_index);
        const auto& stoichiometric_indices = reacDef.getPoolChangeArray(
            Reacdef::PoolChangeArrayType::UPD);
        checkAndUpdateReactionRatesBounds(pMolState, ownerPoint, stoichiometric_indices);
        extent += 1;
    }
    osh::for_each(owned_elements_.data().begin(), owned_elements_.data().end(), [&](osh::LO e) {
        mesh::element_id elem(e);
        for (container::specie_id s(0); s < pMolState.numSpecies(elem); ++s) {
            const auto update_time_interval = period -
                                              pDiffusionVariables.last_update_time(elem, s);
            pDiffusionVariables.occupancy(elem, s) +=
                update_time_interval * static_cast<Omega_h::Real>(pMolState(elem, s));
        }
    });
}

//---------------------------------------------------------

template <typename RNG>
kproc::Reactions::index_type RSSAOperator<RNG>::singleOutReaction() {
    PetscScalar selector = a_upper_bound_.sum() * uniform_(rng_), partial_sums = 0.0;
    for (size_t k = 0; k < a_upper_bound_.size(); k++) {
        PetscScalar rate(a_upper_bound_[k]);
        if (rate == 0.0) {
            continue;
        }
        partial_sums += rate;
        if (selector <= partial_sums) {
            return k;
        }
    }
    throw std::invalid_argument("RSSA engine is unable to select a KProc.");
}

// explicit instantiation definitions
template class RSSAOperator<std::mt19937>;
template RSSAOperator<std::mt19937>::RSSAOperator(MolState& mol_state,
                                                  DiffusionVariables& dv,
                                                  KProcState& kproc_state,
                                                  const OmegaHMesh<2>& mesh,
                                                  std::mt19937& t_rng);
template RSSAOperator<std::mt19937>::RSSAOperator(MolState& mol_state,
                                                  DiffusionVariables& dv,
                                                  KProcState& kproc_state,
                                                  const OmegaHMesh<3>& mesh,
                                                  std::mt19937& t_rng);

//---------------------------------------------------------

}  // namespace zee
