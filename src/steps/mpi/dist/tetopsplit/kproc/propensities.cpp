#include "propensities.hpp"

#include "util/debug.hpp"

namespace steps::dist::kproc {

template <unsigned int Policy>
Propensities<Policy>::Propensities()
    : uniform_(0.0, 1.0 - 2.0 * std::numeric_limits<osh::Real>::epsilon()) {}

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

template <unsigned int Policy>
void Propensities<Policy>::init(const std::array<unsigned, num_kproc_types()>& k_proc_ty) {
    k_proc_ty_2_num_k_proc_ = k_proc_ty;
    // a2ab is a mapping from the kproc type to the last index plus one of
    // propensities of all kprocs related to this kproc type
    std::partial_sum(k_proc_ty_2_num_k_proc_.begin(), k_proc_ty_2_num_k_proc_.end(), a2ab_.begin());
}

template <unsigned int PolicyF>
std::ostream& operator<<(
    std::ostream& ostr,
    const PropensitiesGroup<PolicyF, std::enable_if_t<PropensitiesTraits<PolicyF>::is_direct>>&
        pg) {
    return ostr << "PropensityGroup (Direct)\n"
                << "  idx_: " << pg.idx_ << "\n  partial_sums_: " << pg.partial_sums_;
}

template <unsigned int PolicyF>
std::ostream& operator<<(
    std::ostream& ostr,
    const PropensitiesGroup<PolicyF,
                            std::enable_if_t<PropensitiesTraits<PolicyF>::is_gibson_bruck>>& pg) {
    return ostr << "PropensityGroup (GibsonBruck)" << '\n' << pg.events_;
}

// explicit template instantiation definitions
template class Propensities<PropensitiesPolicy::direct_without_next_event>;
template class Propensities<PropensitiesPolicy::gibson_bruck_without_next_event>;
template class Propensities<PropensitiesPolicy::direct_with_next_event>;
template class Propensities<PropensitiesPolicy::gibson_bruck_with_next_event>;

template struct PropensitiesGroup<PropensitiesPolicy::direct_without_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_without_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::direct_with_next_event>;
template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_with_next_event>;

}  // namespace steps::dist::kproc
