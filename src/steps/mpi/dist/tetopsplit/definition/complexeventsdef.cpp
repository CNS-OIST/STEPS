/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 3,
#    as published by the Free Software Foundation.
#
#    STEPS is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program. If not, see <http://www.gnu.org/licenses/>.
#
#################################################################################   

 */

// STL headers.
#include <algorithm>
#include <boost/smart_ptr/shared_ptr.hpp>
#include <string>

// STEPS headers.
#include "complexeventsdef.hpp"
#include "mpi/dist/tetopsplit/mol_state.hpp"
#include "statedef.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"
#include "util/vocabulary.hpp"
#include "wmdirect/complexevents.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::dist {

////////////////////////////////////////////////////////////////////////////////

ComplexEventdef::ComplexEventdef(const steps::model::ComplexEvent& ce, const Statedef& sd)
    : pcomplexIdx(sd.getComplexModelIdx(model::complex_name(ce.complexId()))) {}

////////////////////////////////////////////////////////////////////////////////

ComplexCreateEventdef::ComplexCreateEventdef(const steps::model::ComplexCreateEvent& ce,
                                             const Statedef& sd)
    : ComplexEventdef::ComplexEventdef(ce, sd)
    , pinit(ce.init()) {}

std::set<model::complex_substate_id> ComplexCreateEventdef::getUpdSet() const {
    std::set<model::complex_substate_id> updset;
    for (auto sus: pinit.range()) {
        if (pinit[sus] > 0) {
            updset.insert(sus);
        }
    }
    return updset;
}

////////////////////////////////////////////////////////////////////////////////

ComplexLHSEventdef::ComplexLHSEventdef(const steps::model::ComplexLHSEvent& ce, const Statedef& sd)
    : ComplexEventdef::ComplexEventdef(ce, sd) {
    pfilters.reserve(ce.filters().size());
    for (const auto& f: ce.filters()) {
        pfilters.emplace_back(f);
    }
}

bool ComplexLHSEventdef::isSame(const std::shared_ptr<const ComplexEventdef>& ev) const {
    const auto* ev2 = dynamic_cast<const ComplexLHSEventdef*>(ev.get());
    return ev2 != nullptr and ComplexEventdef::isSame(ev) and pfilters == ev2->filters();
}

std::set<model::complex_substate_id> ComplexLHSEventdef::getDepSet() const {
    std::set<model::complex_substate_id> depset;
    for (auto& filt: pfilters) {
        for (auto sus: filt.range()) {
            // The event depends on a subunitstate if at least one of its filters depends on it.
            // If the filter does not match the whole possible range of subunitstate counts, it
            // depends on it because the filter could stop matching if the count of the subunitstate
            // drops below filt[sus].min or goes above filt[sus].max.
            if (filt[sus].min > 0 or filt[sus].max < steps::model::COMPLEX_FILTER_MAX_VALUE) {
                depset.insert(sus);
            }
        }
    }
    return depset;
}

////////////////////////////////////////////////////////////////////////////////

ComplexUpdateEventdef::ComplexUpdateEventdef(const steps::model::ComplexUpdateEvent& ce,
                                             const Statedef& sd)
    : ComplexLHSEventdef::ComplexLHSEventdef(ce, sd)
    , preactants(ce.reactants())
    , pupdates(ce.updates())
    , pdestLoc(ce.destLoc())
    , pNoReactants(
          std::all_of(preactants.begin(), preactants.end(), [](auto v) { return v == 0; })) {}

double ComplexUpdateEventdef::rateMult(
    const util::strongid_vector<model::complex_substate_id, uint>& state) const {
    double res = 0;
    for (auto& filt: pfilters) {
        double rmult = 1.0;
        for (auto sus: filt.range()) {
            if (state[sus] < filt[sus].min or state[sus] > filt[sus].max) {
                rmult = 0;
                break;
            }
            uint available = state[sus] - (filt[sus].min - preactants[sus]);
            for (uint n = 0; n < preactants[sus]; ++n) {
                // The division corresponds to dividing by factorial preactants[sus] to correct for
                // all possible orderings of identical reactants
                rmult *= static_cast<double>(available - n) / static_cast<double>(n + 1);
            }
        }
        res = std::max(res, rmult);
    }
    return res;
}

bool ComplexUpdateEventdef::isSame(const std::shared_ptr<const ComplexEventdef>& ev) const {
    const auto* ev2 = dynamic_cast<const ComplexUpdateEventdef*>(ev.get());
    return ev2 != nullptr and ComplexLHSEventdef::isSame(ev) and preactants == ev2->reactants() and
           pupdates == ev2->updates() and pdestLoc == ev2->destLoc();
}

bool ComplexUpdateEventdef::sameReactants(const std::shared_ptr<const ComplexEventdef>& ev) const {
    const auto* ev2 = dynamic_cast<const ComplexUpdateEventdef*>(ev.get());
    return ev2 != nullptr and ComplexEventdef::isSame(ev) and preactants == ev2->reactants();
}

std::set<model::complex_substate_id> ComplexUpdateEventdef::getUpdSet() const {
    std::set<model::complex_substate_id> updset;
    for (auto upd: pupdates) {
        for (auto sus: upd.update.range()) {
            if (upd.update[sus] != 0) {
                updset.emplace(sus.get());
            }
        }
    }
    return updset;
}

const util::strongid_vector<model::complex_substate_id, int>& ComplexUpdateEventdef::getUpdate(
    const util::strongid_vector<model::complex_substate_id, uint>& state,
    rng::RNG& rng) const {
    uint ind = 0;
    if (pupdates.size() > 1) {
        ind = rng.get() % pupdates.size();
    }
    for (uint i = 0; i < pupdates.size(); ++i) {
        const auto& upd = pupdates[(ind + i) % pupdates.size()];
        bool ok = true;
        for (auto sus: upd.requirement.range()) {
            if (upd.requirement[sus] > state[sus]) {
                ok = false;
                break;
            }
        }
        if (ok) {
            return upd.update;
        }
    }
    // Should never reach there since there should always be one element in pupdates that has its
    // requirements satisfied, otherwise the event shouldn't have been selected because it should
    // have 0 rate.
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ComplexDeleteEventdef::ComplexDeleteEventdef(const steps::model::ComplexDeleteEvent& ce,
                                             const Statedef& sd)
    : ComplexLHSEventdef::ComplexLHSEventdef(ce, sd) {}

double ComplexDeleteEventdef::rateMult(
    const util::strongid_vector<model::complex_substate_id, uint>& /*state*/) const {
    return 1.0;
}

std::set<model::complex_substate_id> ComplexDeleteEventdef::getUpdSet() const {
    std::set<model::complex_substate_id> updset;
    for (auto filt: pfilters) {
        for (auto sus: filt.range()) {
            if (filt[sus].max > 0) {
                updset.insert(sus);
            }
        }
        break;
    }
    return updset;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

std::size_t FilterHash::operator()(
    std::vector<util::strongid_vector<model::complex_substate_id,
                                      steps::model::SubunitStateFilter>> const& filts) const {
    std::size_t seed = filts.size();
    for (auto& filt: filts) {
        for (auto& f: filt) {
            seed ^= f.min + 0x9e3779b9 + (seed << 6) + (seed >> 2);
            seed ^= f.max + 0x9e3779b9 + (seed << 6) + (seed >> 2);
        }
    }
    return seed;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

ComplexFilterDescr::ComplexFilterDescr(const steps::model::ComplexFilterDescr& cfd,
                                       const Statedef& sd)
    : complexId(sd.getComplexModelIdx(model::complex_name(cfd.complexId))) {
    filters.reserve(cfd.filters.size());
    for (const auto& f: cfd.filters) {
        filters.emplace_back(f);
    }
}
////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

template <typename Entity>
ComplexLHSCandidates<Entity>::ComplexLHSCandidates(model::complex_id cmplxId, Entity ent)
    : complexIdx_(cmplxId)
    , nbEvents(0)
    , sameReactants(false)
    , filters{nullptr, nullptr}
    , events{nullptr, nullptr}
    , entity(ent) {}

////////////////////////////////////////////////////////////////////////////////

template <typename Entity>
std::vector<std::pair<std::shared_ptr<const ComplexLHSEventdef>, model::complex_individual_id>>
ComplexLHSCandidates<Entity>::selectEvents(const EntityMolecules<Entity>& entmols,
                                           rng::RNG& rng) const {
    std::vector<std::pair<std::shared_ptr<const ComplexLHSEventdef>, model::complex_individual_id>>
        selected;

    double diag = 0;
    double totMult0 = 0, totMult1 = 0;
    double sumCommon0 = 0, sumCommon1 = 0;
    std::map<model::complex_individual_id, double> rateMults0, rateMults1;
    std::unordered_set<model::complex_individual_id> commonCandidates;

    const auto& states = entmols.complexStates(complexIdx_);
    auto [bgn, end] = entmols.complexLocations(complexIdx_).equal_range(entity);
    // For each complex state in the triangle or tetrahedron, check whether it matches the filter of
    // each event. If it matches, compute its rate multiplier and store it, it will be needed later
    // to sample the complex states that are actually selected.
    for (auto it = bgn; it != end; ++it) {
        if (filters[0]->matches(it->second)) {
            auto stit = states.find(it->second);
            if (stit != states.end()) {
                double rm0 = events[0]->rateMult(stit->second);
                rateMults0.emplace(it->second, rm0);
                totMult0 += rm0;
                // If other filter matches, the complex state could be used in either events, but
                // can't be used in both. So we need to keep track of the common candidates.
                if (nbEvents > 1 and filters[1]->matches(it->second)) {
                    double rm1 = events[1]->rateMult(stit->second);
                    rateMults1.emplace(it->second, rm1);
                    commonCandidates.emplace(it->second);
                    totMult1 += rm1;
                    sumCommon0 += rm0;
                    sumCommon1 += rm1;
                    diag += rm0 * rm1;
                }
            }
        } else if (nbEvents > 1 and filters[1]->matches(it->second)) {
            auto stit = states.find(it->second);
            if (stit != states.end()) {
                double rm1 = events[1]->rateMult(stit->second);
                rateMults1.emplace(it->second, rm1);
                totMult1 += rm1;
            }
        }
    }

    // Start with the event that has the least choices
    uint i = 0, j = 1;
    if (nbEvents > 1 and rateMults0.size() > rateMults1.size()) {
        std::swap(i, j);
        std::swap(rateMults0, rateMults1);
        std::swap(totMult0, totMult1);
        std::swap(sumCommon0, sumCommon1);
    }

    if (nbEvents > 1) {
        if (sameReactants) {
            // If the reactants are identical, we need to pay extra attention to the common
            // candidates (the complex states that could be selected for either events). Let us call
            // A and B two complex states that match both events 1 and 2. Let us denote (A:1, B:2)
            // to be the case in which complex state A is matched with event 1 and complex state B
            // is matched with event 2. When reactants for both events are identical, (A:1, B:2) is
            // equivalent to (A:2, B:1). Because of this the total rate multiplier needs to account
            // for this symmetry and avoid counting this cases twice.
            double rmult = totMult0 * totMult1 - (sumCommon0 * sumCommon1 + diag) / 2.0;
            double p = rng.getUnfIE() * rmult;
            double tmp = 0.0;
            bool inCommon = false;
            // Select complex state for the first event (the one with the least choices)
            for (const auto& cand: rateMults0) {
                inCommon = commonCandidates.find(cand.first) != commonCandidates.end();

                if (inCommon) {
                    // If the candidate is common to both event, its probability to be selected for
                    // the first event is proportional to its rate multiplier times the sum of rate
                    // multipliers from all candidates that match the second event but counting only
                    // half of the rate multipliers of the common candidates (to account for the
                    // symmetry mentioned earlier).
                    tmp += cand.second *
                           (totMult1 - (sumCommon1 + rateMults1.at(cand.first)) / 2.0);
                } else {
                    // If the candidate is not common, its probabiliy to be selected is proportional
                    // to its rate multiplier times the sum of rate multipliers from all candidates
                    // that match the second event.
                    tmp += cand.second * totMult1;
                }
                if (tmp >= p) {
                    selected.emplace_back(events[i], cand.first);
                    break;
                }
            }
            // Select second complex state
            if (inCommon) {
                // If the first complex state selected also matched the second event, we need to not
                // select it for the second event.
                model::complex_individual_id ind = selected.back().second;
                p = rng.getUnfIE() * (totMult1 - (sumCommon1 + rateMults1.at(ind)) / 2.0);
                tmp = 0.0;
                for (const auto& cand: rateMults1) {
                    if (cand.first == ind) {
                        // Ignore the complex state that was selected
                        continue;
                    }
                    if (commonCandidates.find(cand.first) != commonCandidates.end()) {
                        // If the complex state is a common candidate, count half its rate
                        // multiplier (because of the symmetry mentioned earlier).
                        tmp += cand.second / 2.0;
                    } else {
                        tmp += cand.second;
                    }
                    if (tmp >= p) {
                        selected.emplace_back(events[j], cand.first);
                        break;
                    }
                }
            } else {
                // If the first complex state selected did not match the second event, we just need
                // to select any candidate that match the second event.
                p = rng.getUnfIE() * totMult1;
                tmp = 0.0;
                for (const auto& cand: rateMults1) {
                    tmp += cand.second;
                    if (tmp >= p) {
                        selected.emplace_back(events[j], cand.first);
                        break;
                    }
                }
            }
        } else {
            // When the reactants are different (A:1, B:2) is not the same as (A:2, B:1) and so we
            // just need to avoid counting cases in which the same complex state would be selected
            // for both events (A:1, A:2)
            double rmult = totMult0 * totMult1 - diag;
            double p = rng.getUnfIE() * rmult;
            double tmp = 0.0;
            bool inCommon = false;
            // Select first complex state
            for (const auto& cand: rateMults0) {
                inCommon = commonCandidates.find(cand.first) != commonCandidates.end();

                if (inCommon) {
                    // If the candidate matches both events, it should be removed from totMult1
                    // (otherwise we would count the case in which it matches both events).
                    tmp += cand.second * (totMult1 - rateMults1.at(cand.first));

                } else {
                    tmp += cand.second * totMult1;
                }
                if (tmp >= p) {
                    selected.emplace_back(events[i], cand.first);
                    break;
                }
            }
            // Select second complex state, avoiding the one that was already selected, as seen
            // above
            if (inCommon) {
                model::complex_individual_id ind = selected.back().second;
                p = rng.getUnfIE() * (totMult1 - rateMults1.at(ind));
                tmp = 0.0;
                for (const auto& cand: rateMults1) {
                    if (cand.first == ind) {
                        continue;
                    }
                    tmp += cand.second;
                    if (tmp >= p) {
                        selected.emplace_back(events[j], cand.first);
                        break;
                    }
                }
            } else {
                p = rng.getUnfIE() * totMult1;
                tmp = 0.0;
                for (const auto& cand: rateMults1) {
                    tmp += cand.second;
                    if (tmp >= p) {
                        selected.emplace_back(events[j], cand.first);
                        break;
                    }
                }
            }
        }
    } else {
        // Only one event
        const double p = rng.getUnfIE() * totMult0;
        double tmp = 0.0;
        for (const auto& cand: rateMults0) {
            tmp += cand.second;
            if (tmp >= p) {
                selected.emplace_back(events[0], cand.first);
                break;
            }
        }
    }

    return selected;
}

template class ComplexLHSCandidates<steps::dist::mesh::tetrahedron_id_t>;
template class ComplexLHSCandidates<steps::dist::mesh::triangle_id_t>;

}  // namespace steps::dist
