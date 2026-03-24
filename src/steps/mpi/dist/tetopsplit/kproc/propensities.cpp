#include "propensities.hpp"

#include "kproc_state.hpp"
#include "mpi/dist/tetopsplit/kproc/kproc_id.hpp"
#include "mpi/dist/tetopsplit/mol_state.hpp"
#include "util/debug.hpp"
#include "util/vocabulary.hpp"

#include <algorithm>
#include <limits>
#include <variant>

namespace steps::dist::kproc {

//-----------------------------------------------

//////////////////
// Propensities //
//////////////////

//-----------------------------------------------

template <unsigned int Policy>
Propensities<Policy>::Propensities()
    : uniform_(0.0, 1.0 - 2.0 * std::numeric_limits<osh::Real>::epsilon()) {}

//-----------------------------------------------

template <unsigned int Policy>
void Propensities<Policy>::init(const std::array<unsigned, num_kproc_types()>& k_proc_ty,
                                const typename propensity_function_traits::value& propensities_fun,
                                const kproc_groups_t& groups) {
    init(k_proc_ty);
    propensities_groups_.clear();
    v_.clear();
    // build a map of kproc type to number of kprocs of type as a vector
    fun_ = propensities_fun;
    // v_ is a placeholder for propensities

    const auto num_propensitites = a2ab_.back();
    v_.resize(num_propensitites, std::numeric_limits<osh::Real>::quiet_NaN());
    propensities_groups_.reserve(static_cast<size_t>(groups.size()));
    local_indices_.resize(num_propensitites);
    for (const auto& group: groups) {
        size_t l{};
        for (auto kp: group) {
            local_indices_[ab(KProcID(static_cast<unsigned>(kp)))] = l;
            l++;
        }
        propensities_groups_.emplace_back(*this, group);
    }
}

//-----------------------------------------------

template <unsigned int Policy>
void Propensities<Policy>::init(const std::array<unsigned, num_kproc_types()>& k_proc_ty) {
    k_proc_ty_2_num_k_proc_ = k_proc_ty;
    // a2ab is a mapping from the kproc type to the last index plus one of
    // propensities of all kprocs related to this kproc type
    std::partial_sum(k_proc_ty_2_num_k_proc_.begin(), k_proc_ty_2_num_k_proc_.end(), a2ab_.begin());
}

//-----------------------------------------------

/////////////////////////////////////
// PropensitiesGroup: Gibson Bruck //
/////////////////////////////////////

//-----------------------------------------------

template <unsigned int Policy>
void PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_gibson_bruck>>::
    update_outdated(const std::vector<KProcID>& outdated,
                    const MolState& mol_state,
                    rng::RNG& rng,
                    const osh::Real state_time) {
    for (auto kid: outdated) {
        adjust_existing_events(kid, mol_state, rng, state_time);
    }
}

//-----------------------------------------------

/////////////////////////////////
// PropensitiesGroup: RLeaping //
/////////////////////////////////

//-----------------------------------------------

template <unsigned int Policy>
void PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_rleaping>>::init(
    const MolState& /*state*/,
    const kproc::KProcState& kpState) {
    // Build groups
    for (size_t i = 0; i < idx_.size(); ++i) {
        updateCRGroups(i);
    }

    std::map<MolStateElementID, size_t> elem2InfoIdx;
    std::vector<std::vector<ReacUpdateInfo>> tmpReacUpdates;
    {
        // Process LHS of KProcs and populate elem2InfoIdx
        auto processKProcLHS = [&, this](const KProcID& kid, const auto& reactions) {
            const auto& lhs = reactions.getMolStateElementsRequiredQuantities(kid.id());
            const auto& rdef = reactions.getReacDef(kid.id());
            unsigned int order = rdef.getOrder();
            for (const auto& elemSpec: lhs) {
                MolStateElementID elemId(std::get<0>(elemSpec), std::get<1>(elemSpec));
                unsigned int mult = -std::get<2>(elemSpec);

                auto [it, added] = elem2InfoIdx.emplace(elemId, specInfos.size());
                if (added) {
                    specInfos.emplace_back(elemId, 0, 0);
                    tmpReacUpdates.emplace_back();
                }
                auto& info = specInfos[it->second];
                if (order > info.h) {
                    info.h = order;
                    info.n = mult;
                } else if (order == info.h and mult > info.n) {
                    info.n = mult;
                }
            }
        };
        for (auto it = ids_.begin(); it != ids_.end(); ++it) {
            KProcID kid(*it);
            switch (kid.type()) {
            case KProcType::Reac:
                processKProcLHS(kid, kpState.reactions());
                break;
            case KProcType::SReac:
                processKProcLHS(kid, kpState.surfaceReactions());
                break;
            case KProcType::VDepSReac:
                processKProcLHS(kid, kpState.vDepSurfaceReactions());
                break;
            case KProcType::GHKSReac:
                processKProcLHS(kid, kpState.ghkSurfaceReactions());
                break;
            default:
                assert(false);
                break;
            }
        }
    }
    {
        // Process updates of KProcs and populate tmpReacUpdates
        auto processKProcUPD =
            [&, this](const KProcID& kid, const size_t& lidx, const auto& reactions) {
                const auto& updates = reactions.getMolStateElementsUpdatesQuantities(kid.id());
                for (const auto& elemSpec: updates) {
                    MolStateElementID elemId(std::get<0>(elemSpec), std::get<1>(elemSpec));
                    auto it = elem2InfoIdx.find(elemId);
                    if (it != elem2InfoIdx.end()) {
                        // Only add the update if the (elem, spec) pair can affect propensities of
                        // kprocs
                        tmpReacUpdates[it->second].emplace_back(lidx, std::get<2>(elemSpec));
                    }
                }
            };
        size_t lidx{};
        for (auto it = ids_.begin(); it != ids_.end(); ++it, ++lidx) {
            KProcID kid(*it);
            switch (kid.type()) {
            case KProcType::Reac:
                processKProcUPD(kid, lidx, kpState.reactions());
                break;
            case KProcType::SReac:
                processKProcUPD(kid, lidx, kpState.surfaceReactions());
                break;
            case KProcType::VDepSReac:
                processKProcUPD(kid, lidx, kpState.vDepSurfaceReactions());
                break;
            case KProcType::GHKSReac:
                processKProcUPD(kid, lidx, kpState.ghkSurfaceReactions());
                break;
            default:
                assert(false);
                break;
            }
        }
    }

    // Fill reacInfos with the computed temporary values
    std::vector<size_t> reacs_per_elemSpec(tmpReacUpdates.size());
    std::transform(tmpReacUpdates.begin(),
                   tmpReacUpdates.end(),
                   reacs_per_elemSpec.begin(),
                   [](const auto& v) { return v.size(); });
    reacInfos.reshape(reacs_per_elemSpec);
    for (unsigned int i = 0; i < tmpReacUpdates.size(); ++i) {
        for (unsigned int j = 0; j < reacs_per_elemSpec[i]; ++j) {
            reacInfos(i, j) = tmpReacUpdates[i][j];
        }
    }
}

//-----------------------------------------------

template <unsigned int Policy>
std::vector<std::pair<KProcID, unsigned int>>
PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_rleaping>>::drawEvents(
    rng::RNG& rng,
    unsigned int L) const {
    std::vector<std::pair<KProcID, unsigned int>> events;
    int sel = 0;
    int indsel = posGroups.size();
    double sum_ai = 0;
    while (L > 0) {
        // Walk through the groups, from highest summed propensity to lowest
        if (not getNextGroup(sel, indsel)) {
            break;
        }
        const auto& group = sel >= 0 ? posGroups[indsel] : negGroups[indsel];
        if (group.indices.empty()) {
            continue;
        }
        osh::Real avgProp = group.sum / static_cast<osh::Real>(group.indices.size());
        if (avgProp / (a0_ - sum_ai) * L >= binomial_threshold()) {
            // Use binomial sampling in the group
            for (unsigned int i = 0; i < group.indices.size() and L > 0; i++) {
                osh::Real ai = propensities_.v_[idx_[group.indices[i]]];
                std::binomial_distribution<> dist(L, ai / (a0_ - sum_ai));
                unsigned int nb = dist(rng);
                if (nb > 0) {
                    events.emplace_back(ids_(group.indices[i]), nb);
                    L -= nb;
                }
                sum_ai += ai;
            }
        } else {
            // Sample all future reactions one by one
            for (; L > 0; --L) {
                osh::Real selector = propensities_.uniform_(rng) * (a0_ - sum_ai);
                int sel2 = sel;
                int indsel2 = indsel;
                osh::Real partial_sum = 0;
#ifndef NDEBUG
                bool found = false;
#endif
                do {
                    const auto& group2 = sel2 >= 0 ? posGroups[indsel2] : negGroups[indsel2];
                    if (group2.indices.empty()) {
                        continue;
                    }
                    partial_sum += group2.sum;
                    // There might be a potential issue here since we are summing group
                    // propensities from highest to smallest. It could be that the selector is
                    // always bigger than the partial_sum.
                    // Although this possibility is dealt with in the serial implementation of
                    // composition-rejection, it is unlikely enough in the distributed case that
                    // we do not tackle it here.
                    if (partial_sum >= selector) {
                        // Do rejection sampling in that group
                        int rnd_ind = rng.get() % group2.indices.size();
                        osh::Real val = rng.getUnfII() *
                                        (sel2 >= 0 ? (1 << indsel2) : (1 >> -indsel2));
                        while (propensities_.v_[idx_[group2.indices[rnd_ind]]] < val) {
                            rnd_ind = rng.get() % group2.indices.size();
                            val = rng.getUnfII() * (sel2 >= 0 ? (1 << indsel2) : (1 >> -indsel2));
                        }
                        events.emplace_back(ids_(group2.indices[rnd_ind]), 1);
#ifndef NDEBUG
                        found = true;
#endif
                        break;
                    }
                } while (getNextGroup(sel2, indsel2));
                // Check that a reaction could be found
                assert(found);
            }
        }
    }
    assert(L == 0);
    return events;
}

//-----------------------------------------------

template <unsigned int Policy>
long long
PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_rleaping>>::computeL(
    const MolState& state,
    const double& epsilon,
    const double& theta) const {
    if (a0_ <= 0) {
        return 0;
    }
    double inva0 = 1.0 / a0_;
    double value = std::numeric_limits<double>::infinity();
    for (unsigned int i = 0; i < specInfos.size(); ++i) {
        const auto& info = specInfos[i];
        double xi = state(info.elemId);
        // Compute g_i
        double gi = info.h;
        double sum = 0;
        for (double k = 1; k < info.n; k += 1.0) {
            sum += k / (xi - k);
        }
        gi += static_cast<double>(info.h) / static_cast<double>(info.n) * sum;
        // Compute \mu_i and \sigma_i
        double mu = 0;
        double sigma2 = 0;
        for (const auto& reacInfo: reacInfos[i]) {
            const auto& ai = propensities_.v_[idx_[reacInfo.ridx]];
            mu += static_cast<double>(reacInfo.update) * ai;
            sigma2 += static_cast<double>(reacInfo.update * reacInfo.update) * ai;
            // Non-negative species control
            if (reacInfo.update < 0 and ai > std::numeric_limits<double>::epsilon()) {
                // Do the computation with everything divided by a0
                value = std::min(value,
                                 (inva0 - theta * (inva0 - 1.0 / ai)) *
                                     static_cast<double>(static_cast<int>(xi) / -reacInfo.update));
            }
        }
        double num = std::max(epsilon * xi / gi, 1.0);
        value = std::min(value,
                         std::min(num / std::abs(mu),
                                  num * num / (std::abs(sigma2) - std::abs(mu * mu * inva0))));
        if (value < inva0) {
            return 1;
        }
    }
    return std::max(static_cast<long long>(std::floor(a0_ * value)), 1LL);
}

//-----------------------------------------------

template <unsigned int PolicyF>
std::ostream& operator<<(
    std::ostream& ostr,
    const PropensitiesGroup<PolicyF, std::enable_if_t<PropensitiesTraits<PolicyF>::is_direct>>&
        pg) {
    return ostr << "PropensityGroup (Direct)\n"
                << "  idx_: " << pg.idx_ << "\n  partial_sums_: " << pg.partial_sums_;
}

//-----------------------------------------------

template <unsigned int PolicyF>
std::ostream& operator<<(
    std::ostream& ostr,
    const PropensitiesGroup<PolicyF,
                            std::enable_if_t<PropensitiesTraits<PolicyF>::is_gibson_bruck>>& pg) {
    return ostr << "PropensityGroup (GibsonBruck)" << '\n' << pg.events_;
}

//-----------------------------------------------

// explicit template instantiation definitions
template class Propensities<PropensitiesPolicy::direct_without_next_event>;
template class Propensities<PropensitiesPolicy::gibson_bruck_without_next_event>;
template class Propensities<PropensitiesPolicy::direct_with_next_event>;
template class Propensities<PropensitiesPolicy::gibson_bruck_with_next_event>;
template class Propensities<PropensitiesPolicy::rleaping_without_next_event>;

template struct PropensitiesGroup<PropensitiesPolicy::direct_without_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_without_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::direct_with_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_with_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::rleaping_without_next_event>;

}  // namespace steps::dist::kproc
