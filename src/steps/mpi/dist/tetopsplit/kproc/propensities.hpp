#pragma once

#include <array>
#include <iosfwd>
#include <random>
#include <vector>

#include "event_queue.hpp"
#include "kproc_id.hpp"
#include "mpi/dist/tetopsplit/fwd.hpp"
#include "mpi/dist/tetopsplit/mol_state.hpp"
#include "rng/rng.hpp"
#include "util/error.hpp"
#include "util/flat_multimap.hpp"
#include "util/vocabulary.hpp"

namespace steps::dist::kproc {

/**
 * Bitmask values to describe at compile-time the behavior of \a Propensities
 * and \a PropensitiesGroup using template class specialization based on the value of the
 * bitmask.
 */
struct PropensitiesPolicy {
    static constexpr unsigned int without_next_event = 0b1;
    static constexpr unsigned int with_next_event = 0b10;
    static constexpr unsigned int direct_event = 0b100;
    static constexpr unsigned int gibson_bruck_event = 0b1000;

    static constexpr unsigned int direct_without_next_event = without_next_event | direct_event;
    static constexpr unsigned int direct_with_next_event = with_next_event | direct_event;
    static constexpr unsigned int gibson_bruck_without_next_event = without_next_event |
                                                                    gibson_bruck_event;
    static constexpr unsigned int gibson_bruck_with_next_event = with_next_event |
                                                                 gibson_bruck_event;

    static constexpr unsigned int default_policy = direct_with_next_event;

    static constexpr unsigned int search_method_mask = direct_event | gibson_bruck_event;
    static constexpr unsigned int next_event_mask = with_next_event | without_next_event;

    /**
     * constexpr function (evaluated at compile-time) to translate a \a NextEventSearchMethod
     * enumerated value to its corresponding bitmask.
     *
     * This is the C++17-way to write a type traits. The previous way could have been:
     * \code
     * template <NextEventSearchMethod SearchMethod> class GetPropensitiesPolicy {
     *     static_assert(false, "Missing specialization for this method");
     * };
     * template <>
     * class GetPropensitiesPolicy<NextEventSearchMethod::Direct> {
     *     static constexpr bool value = PropensitiesPolicy::direct_event;
     * };
     * template <>
     * class GetPropensitiesPolicy<NextEventSearchMethod::GibsonBruck> {
     *     static constexpr bool value = PropensitiesPolicy::gibson_bruck_event;
     * };
     * \endcode
     */
    template <NextEventSearchMethod SearchMethod>
    static constexpr unsigned int get() {
        switch (SearchMethod) {
        case NextEventSearchMethod::Direct:
            return direct_event;
        case NextEventSearchMethod::GibsonBruck:
            return gibson_bruck_event;
        default:
            static_assert(true, "Unexpected enum value");
        }
    }
};

/**
 *
 * \tparam Policy bitmask to extract information from
 */
template <unsigned int Policy>
struct PropensitiesTraits {
    static_assert((Policy & PropensitiesPolicy::search_method_mask) != 0,
                  "a search method must be specified");
    static_assert((Policy & PropensitiesPolicy::search_method_mask) !=
                      PropensitiesPolicy::search_method_mask,
                  "only one search method must be specified");
    static_assert((Policy & PropensitiesPolicy::next_event_mask) != 0,
                  "event management must be specified");
    static_assert((Policy & PropensitiesPolicy::next_event_mask) !=
                      PropensitiesPolicy::next_event_mask,
                  "only one event manager must be specified");

    /// true if the Gibson-Bruck search method is selected, false otherwise
    static constexpr bool is_gibson_bruck = (Policy & PropensitiesPolicy::gibson_bruck_event) != 0;
    /// true if the direect method is selected, false otherwise
    static constexpr bool is_direct = (Policy & PropensitiesPolicy::direct_event) != 0;
    /// true if the propensities should handle the next event, false otherwise
    static constexpr bool handle_next_event = (Policy & PropensitiesPolicy::with_next_event) != 0;
};

template <unsigned Policy = PropensitiesPolicy::default_policy, class Enable = void>
struct PropensitiesGroup;

using EventTime = osh::Real;
using Event = std::pair<EventTime, kproc::KProcID>;

/// Hold all indep. kproc
using kproc_groups_t = util::flat_multimap<osh::LO, 1>;
using kproc_group_t = kproc_groups_t::const_element_type;

/// Hold all dependenc
/// ies of a given kinetic process
using dependencies_t = util::flat_multimap<osh::LO, 1>;
/// the kinetic processes that depend on a change of propensity
using KProcDeps = dependencies_t::const_element_type;

/**
 * \brief Manage propensities of all kprocs into a flat multimap structure.
 * \tparam Policy mask
 *   * propgpr_with_next_event: if true, propensities groups can generate a next event.
 */
template <unsigned int Policy = PropensitiesPolicy::default_policy>
class Propensities {
  public:
    friend struct PropensitiesGroup<Policy>;
    friend class KProcState;

    static_assert(
        std::is_same<unsigned, typename std::underlying_type<kproc::KProcType>::type>::value,
        "KProcType needs to be unsigned.");

    Propensities();

    /**
     * \brief Initialize the propensities of classes of kprocs.
     *
     * \param k_proc_ty a vector of kprocs types handled alongside number of such kprocs
     * \param propensities_fun a function to compute the propensity of a kproc
     */
    void init(const std::array<unsigned, kproc::num_kproc_types()>& k_proc_ty,
              const typename propensity_function_traits::value& propensities_fun,
              const kproc_groups_t& groups);

    /**
     * \brief Query the propensity of a given index.
     *
     * \param idx an index
     * \return the propensity of the kproc idx
     */
    inline osh::Real operator[](size_t idx) const noexcept {
        return v_[idx];
    }

    /**
     * \brief Query the propensity of a given kproc.
     *
     * \param id kroc id
     * \return the propensity of the kproc
     */
    inline osh::Real operator[](KProcID id) const noexcept {
        return v_[ab(id)];
    }

    /**
     * \brief Query the kproc id associated with a propensity index.
     *
     * \param idx propensity index
     * \return the kproc id
     */
    KProcID kProcId(size_t idx) const {
        auto idx_l = std::distance(a2ab_.begin(),
                                   std::upper_bound(a2ab_.begin(), a2ab_.end(), idx));
        assert(idx_l >= 0);
        auto idx_ul = static_cast<unsigned>(idx_l);
        assert(idx_ul < a2ab_.size());
        if (idx_ul > 0) {
            return {KProcType(idx_ul), static_cast<unsigned>(idx - a2ab_[idx_ul - 1])};
        } else {
            return {KProcType(0), static_cast<unsigned>(idx)};
        }
    }

    /**
     *
     * \return the number of all kprocs
     */
    inline size_t size() const noexcept {
        return a2ab_.back();
    }

    /**
     * \return all propensities groups
     */
    const auto& groups() const noexcept {
        return propensities_groups_;
    }

    /**
     * \return all propensities groups
     */
    auto& groups() noexcept {
        return propensities_groups_;
    }

  private:
    /**
     * \brief Setup the flat multimap structure of propensities across the various
     * kinetic processes.
     *
     * \param k_proc_ty a vector of kprocs types handled alongside number of such
     * kprocs
     */
    void init(const std::array<unsigned, num_kproc_types()>& k_proc_ty);

    /**
     * \brief Propensity index of a given kproc.
     *
     * \param kp a proc id
     * \return the index of the propensity
     */
    inline size_t ab(KProcID kp) const noexcept {
        auto kpt = static_cast<size_t>(kp.type());
        assert(kpt < a2ab_.size());
        assert(kp.id() < a2ab_[kpt]);
        if (kpt > 0) {
            return a2ab_[static_cast<size_t>(kp.type()) - 1] + kp.id();
        } else {
            return kp.id();
        }
    }

    std::vector<osh::Real> v_;
    std::vector<size_t> local_indices_;
    std::array<unsigned, num_kproc_types()> a2ab_{};
    std::array<unsigned, num_kproc_types()> k_proc_ty_2_num_k_proc_{};
    typename propensity_function_traits::value fun_;
    std::uniform_real_distribution<double> uniform_;
    std::vector<PropensitiesGroup<Policy>> propensities_groups_;
};

//--------------------------------------------------------

/**
 * \brief A group of propensities where next event is searched via Gibson Bruck
 * algorithm. Gibson, B. and Bruck, J. “Efficient Exact Stochastic Simulation of
 * Chemical Systems with Many Species and Many Channels.” The Journal of
 * Physical Chemistry A, 104 9 (2000): 1876-1889.
 *
 */
template <unsigned int Policy>
struct PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_gibson_bruck>> {
    /**
     * \brief Ctor.
     *
     * \param propensities propensities
     * \param ids KProcIds of processes handled by the group
     */
    PropensitiesGroup(Propensities<Policy>& propensities, kproc_group_t ids)
        : kprocs_(ids)
        , propensities_(propensities) {}

    static constexpr bool handle_next_event() {
        return PropensitiesTraits<Policy>::handle_next_event;
    }

    inline osh::Real getExp(double prop, rng::RNG& rng) const {
        return rng.getExp(prop);
    }

    /**
     * \brief reset the group data structure
     */
    void reset(const MolState& mol_state, rng::RNG& rng, const osh::Real state_time) {
        events_.clear();
        for (auto kp: kprocs_) {
            KProcID kid{static_cast<unsigned>(kp)};
            auto idx = propensities_.ab(kid);
            osh::Real propensity = propensities_.fun_(kid, mol_state);
            propensities_.v_[idx] = propensity;
            if constexpr (!handle_next_event()) {
                continue;
            }
            if (propensity > std::numeric_limits<osh::Real>::epsilon()) {
                events_.update(kid, state_time + getExp(propensity, rng));
            } else {
                events_.update(kid, std::numeric_limits<osh::Real>::infinity());
            }
        }
    }

    /**
     * \brief Update the maximum time threshold of the priority queue
     * in the Gibson Bruck method.
     *
     * An event is added to the next event queue only if it
     * happens no later than max_time, otherwise it is stored in the reserve.
     * When this function is called, the solution loops over all events stored
     * in the reserve and inserts the ones below or equal to the new threshold
     * to the next event queue.
     *
     * \param max_time the new maximum time threshold
     */
    void updateMaxTime(const osh::Real max_time) {
        events_.updateMaxTime(max_time);
    }

    /**
     * \brief Update the state of all the propensities and next_event structure following Gibson
     * Bruck algorithm.
     *
     * \param mol_state molecular state
     * \param rng random number generator
     */
    void update_all(MolState& mol_state, rng::RNG& rng, const osh::Real state_time) {
        for (auto kp: kprocs_) {
            KProcID kid{static_cast<unsigned>(kp)};
            adjust_existing_events(kid, mol_state, rng, state_time);
        }
    }

    /**
     * \brief Update the state of a selected number of propensities (only tets with diffusing
     * molecules need to update propensities) and next_event structure following Gibson Bruck
     * algorithm.
     *
     * \param mol_state molecular state
     * \param rng random number generator
     */
    void update_outdated(MolState& mol_state, rng::RNG& rng, const osh::Real state_time) {
        for (auto kp: mol_state.outdated_kprocs()) {
            KProcID kid{kp};
            adjust_existing_events(kid, mol_state, rng, state_time);
        }
        mol_state.outdated_kprocs().clear();
    }

    /**
     * \brief Update Gibson Bruck next reaction structure following the
     * occurrence of a kinetic event.
     *
     * \param mol_state molecular state
     * \param rng random generator
     * \param event a kinetic event that occurred
     * \param selection a selection of kprocs that need recomputation following
     * the event
     */
    template <typename T>
    void update(const MolState& mol_state, rng::RNG& rng, const Event event, const T& selection) {
        using cast_type =
            typename std::conditional<std::is_same<T, KProcDeps>::value, unsigned, KProcID>::type;
        if constexpr (!handle_next_event()) {
            for (auto k: selection) {
                KProcID kp(static_cast<cast_type>(k));
                auto idx = propensities_.ab(kp);
                auto new_propensity = propensities_.fun_(kp, mol_state);
                propensities_.v_[idx] = new_propensity;
            }
        } else {
            // process the event
            {
                const auto& kp = event.second;
                auto idx = propensities_.ab(kp);
                const osh::Real new_propensity = propensities_.fun_(kp, mol_state);
                propensities_.v_[idx] = new_propensity;
                if (new_propensity > std::numeric_limits<osh::Real>::epsilon()) {
                    events_.update(kp, event.first + getExp(new_propensity, rng));
                } else {
                    events_.update(kp, std::numeric_limits<osh::Real>::infinity());
                }
            }

            // process its dependencies
            for (auto k: selection) {
                KProcID kp(static_cast<cast_type>(k));
                if (kp.data() != event.second.data()) {
                    adjust_existing_events(kp, mol_state, rng, event.first);
                }
            }
        }
    }

    void adjust_existing_events(KProcID kp,
                                const MolState& mol_state,
                                rng::RNG& rng,
                                const osh::Real current_state_time) {
        auto idx = propensities_.ab(kp);
        // unrelated to the event that's just happened
        osh::Real& old_propensity = propensities_.v_[idx];
        osh::Real new_propensity = propensities_.fun_(kp, mol_state);
        if (old_propensity != new_propensity) {
            if (new_propensity > std::numeric_limits<osh::Real>::epsilon()) {
                osh::Real adj_time;
                auto old_time = events_.getEventTime(kp);
                if (old_propensity > std::numeric_limits<osh::Real>::epsilon() &&
                    old_time != current_state_time) {
                    assert(old_time > current_state_time);
                    adj_time = old_propensity / new_propensity * (old_time - current_state_time) +
                               current_state_time;
                } else {
                    adj_time = current_state_time + getExp(new_propensity, rng);
                }
                old_propensity = new_propensity;
                events_.update(kp, adj_time);
            } else {
                old_propensity = new_propensity;
                events_.update(kp, std::numeric_limits<osh::Real>::infinity());
            }
        }
    }
    /**
     * \brief Draw a kproc id from a discrete distribution of probabilities
     * given by scaled propensities.
     *
     * \return a kproc sample
     */
    Event drawEvent(rng::RNG& /*rng*/, osh::Real /* sim_time */) {
        return events_.getFirst();
    }

    ///  pretty printer
    template <unsigned int PolicyF>
    friend std::ostream& operator<<(
        std::ostream& os,
        const PropensitiesGroup<PolicyF,
                                std::enable_if_t<PropensitiesTraits<PolicyF>::is_gibson_bruck>>&
            pg);

  private:
    EventQueue events_;
    kproc_group_t kprocs_;
    Propensities<Policy>& propensities_;
};

//--------------------------------------------------------

/**
 * \brief A group of propensities where next event is searched via the Direct
 * method of Gillespie. D. Gillespie, Stochastic Simulation of Chemical
 * Kinetics, Annu. Rev. Phys. Chem. 2007. 58:35-55
 *
 */
template <unsigned int Policy>
struct PropensitiesGroup<Policy, std::enable_if_t<PropensitiesTraits<Policy>::is_direct>> {
    /**
     * \brief Ctor.
     *
     * \param propensities all propensities of kprocs
     * \param ids KProcIds of kprocs handled by the group
     */
    PropensitiesGroup(Propensities<Policy>& propensities, const kproc_group_t& ids)
        : idx_(static_cast<size_t>(ids.size()))
        , ids_(ids)
        , propensities_(propensities) {
        std::transform(ids.begin(), ids.end(), idx_.begin(), [&propensities](osh::LO id) {
            return propensities.ab(KProcID(static_cast<unsigned>(id)));
        });
        if constexpr (handle_next_event()) {
            partial_sums_.resize(static_cast<size_t>(ids.size()),
                                 std::numeric_limits<osh::Real>::quiet_NaN());
        }
        idx_.shrink_to_fit();
        partial_sums_.shrink_to_fit();
    }

    static constexpr bool handle_next_event() {
        return PropensitiesTraits<Policy>::handle_next_event;
    }

    /**
     * \brief reset the group data structure
     */
    void reset(const MolState& /*mol_state*/, rng::RNG& /*rng*/, const osh::Real /*state_time*/) {
        // do nothing
    }

    void updateMaxTime(const osh::Real /*max_time*/) {
        // do nothing
    }

    /**
     * \brief Update the propensities of all kprocs and state of the class.
     *
     * \param mol_state molecular state
     */
    void update_all(const MolState& mol_state, rng::RNG& /*rng*/, const osh::Real /*state_time*/) {
        size_t k{};
        for (auto it = ids_.begin(); it != ids_.end(); it++, k++) {
            propensities_.v_[idx_[k]] = propensities_.fun_(KProcID(static_cast<unsigned>(*it)),
                                                           mol_state);
        }
        if constexpr (handle_next_event()) {
            updatePartialSums(0);
        }
    }

    /**
     * \brief Update only the outdated propensities and state of the class.
     *
     * \param mol_state molecular state
     */
    void update_outdated(const MolState& mol_state,
                         rng::RNG& /*rng*/,
                         const osh::Real /*state_time*/) {
        /* TODO : Follow the same approach as in Gibson Bruck, i.e. update outdated propensities
           only (test needed). For now go with the original version, i.e. update all propensities.
         */
        size_t k{};
        for (auto it = ids_.begin(); it != ids_.end(); it++, k++) {
            propensities_.v_[idx_[k]] = propensities_.fun_(KProcID(static_cast<unsigned>(*it)),
                                                           mol_state);
        }
        if constexpr (handle_next_event()) {
            updatePartialSums(0);
        }
    }

    /**
     * \brief Update the state of propensities of a selected number of kprocs.
     *
     * \param mol_state molecular state
     * \param selection of kprocs that need update in the current group
     */
    template <typename T>
    void update(const MolState& mol_state,
                rng::RNG& /*rng*/,
                const Event& /*event*/,
                const T& selection) {
        using cast_type =
            typename std::conditional<std::is_same<T, KProcDeps>::value, unsigned, KProcID>::type;
        size_t min_idx(partial_sums_.size());
        for (auto k: selection) {
            KProcID kp(static_cast<cast_type>(k));
            min_idx = std::min(min_idx, propensities_.local_indices_[propensities_.ab(kp)]);
            auto idx = propensities_.ab(kp);
            propensities_.v_[idx] = propensities_.fun_(kp, mol_state);
        }
        if constexpr (handle_next_event()) {
            updatePartialSums(min_idx);
        }
    }

    /**
     * \brief Update propensities partial sums in the given group
     */
    void updatePartialSums(size_t min_idx) noexcept {
        osh::Real val{};
        if (min_idx > 0) {
            val = partial_sums_[min_idx - 1];
        }
        for (size_t k = min_idx; k < partial_sums_.size(); k++) {
            val += propensities_[idx_[k]];
            partial_sums_[k] = val;
        }
    }

    /**
     * \brief Draw a kproc id from a discrete distribution of probabilities
     * given by scaled propensities.
     *
     * There are 2 somewhat independent drawing here:
     *
     * - when the next event occurs:
     *    it is a roll on an exponential distribution with a characteristic value equal to
     *    the sum of the propensities of the group
     * - what proc is the next event:
     *    it is a uniform roll over the partial sums, The bigger your weight, the higher is the
     * probability of being picked
     *
     * \param rng a random number generator
     * \return a kproc sample
     */
    Event drawEvent(rng::RNG& rng, osh::Real sim_time) {
        // when the propensities are 0 or there we assume that the next event is at infinity
        if (partial_sums_.empty()) {
            return {std::numeric_limits<osh::Real>::infinity(), KProcID(0)};
        }
        if (partial_sums_.back() < std::numeric_limits<osh::Real>::epsilon()) {
            return {std::numeric_limits<osh::Real>::infinity(), KProcID(0)};
        }
        const auto next_arrival = rng.getExp(partial_sums_.back());
        auto idx =
            std::distance(partial_sums_.begin(),
                          std::upper_bound(partial_sums_.begin(),
                                           partial_sums_.end(),
                                           partial_sums_.back() * propensities_.uniform_(rng)));
        assert(idx >= 0 && static_cast<size_t>(idx) < idx_.size());
        return {sim_time + next_arrival, propensities_.kProcId(idx_[static_cast<size_t>(idx)])};
    }

    /// pretty printer
    template <unsigned int PolicyF>
    friend std::ostream& operator<<(
        std::ostream& ostr,
        const PropensitiesGroup<PolicyF, std::enable_if_t<PropensitiesTraits<PolicyF>::is_direct>>&
            pg);

  private:
    std::vector<size_t> idx_;
    kproc::kproc_group_t ids_;
    std::vector<osh::Real> partial_sums_;
    Propensities<Policy>& propensities_;
};

/**
 * \a PropensitiesGroup pretty printer
 */
template <unsigned int Policy>
inline std::ostream& operator<<(std::ostream& ostr, const PropensitiesGroup<Policy>& /* group */) {
    return ostr << "PropensityGroup (base)\n";
}

//--------------------------------------------------------

// explicit template instantiation declarations
extern template class Propensities<PropensitiesPolicy::direct_without_next_event>;
extern template class Propensities<PropensitiesPolicy::gibson_bruck_without_next_event>;
extern template class Propensities<PropensitiesPolicy::direct_with_next_event>;
extern template class Propensities<PropensitiesPolicy::gibson_bruck_with_next_event>;

extern template struct PropensitiesGroup<PropensitiesPolicy::direct_without_next_event>;
extern template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_without_next_event>;
extern template struct PropensitiesGroup<PropensitiesPolicy::direct_with_next_event>;
extern template struct PropensitiesGroup<PropensitiesPolicy::gibson_bruck_with_next_event>;

}  // namespace steps::dist::kproc
