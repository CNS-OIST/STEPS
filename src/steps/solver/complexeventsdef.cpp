/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2023 Okinawa Institute of Science and Technology, Japan.
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

#include "solver/complexeventsdef.hpp"

#include "rng/rng.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"

namespace steps::solver {

ComplexEventdef::ComplexEventdef(const model::ComplexEvent& ce, const Statedef& sd)
    : pcomplexIdx(sd.getComplexIdx(ce.complexId())) {}

////////////////////////////////////////////////////////////////////////////////

ComplexCreateEventdef::ComplexCreateEventdef(const model::ComplexCreateEvent& ce,
                                             const Statedef& sd)
    : ComplexEventdef::ComplexEventdef(ce, sd)
    , pinit(ce.init()) {}

std::set<complex_substate_id> ComplexCreateEventdef::getUpdSet() const {
    std::set<complex_substate_id> updset;
    for (auto sus: pinit.range()) {
        if (pinit[sus] > 0) {
            updset.insert(sus);
        }
    }
    return updset;
}

////////////////////////////////////////////////////////////////////////////////

ComplexLHSEventdef::ComplexLHSEventdef(const model::ComplexLHSEvent& ce, const Statedef& sd)
    : ComplexEventdef::ComplexEventdef(ce, sd)
    , pfilters(ce.filters()) {}

bool ComplexLHSEventdef::isSame(const std::shared_ptr<const ComplexEventdef>& ev) const {
    const auto* ev2 = dynamic_cast<const ComplexLHSEventdef*>(ev.get());
    return ev2 != nullptr and ComplexEventdef::isSame(ev) and pfilters == ev2->filters();
}

std::set<complex_substate_id> ComplexLHSEventdef::getDepSet() const {
    std::set<complex_substate_id> depset;
    for (auto& filt: pfilters) {
        for (auto sus: filt.range()) {
            if (filt[sus].min > 0 or filt[sus].max < model::COMPLEX_FILTER_MAX_VALUE) {
                depset.insert(sus);
            }
        }
    }
    return depset;
}

////////////////////////////////////////////////////////////////////////////////

ComplexUpdateEventdef::ComplexUpdateEventdef(const model::ComplexUpdateEvent& ce,
                                             const Statedef& sd)
    : ComplexLHSEventdef::ComplexLHSEventdef(ce, sd)
    , preactants(ce.reactants())
    , pupdates(ce.updates())
    , pdestLoc(ce.destLoc()) {}

double ComplexUpdateEventdef::rateMult(
    const util::strongid_vector<complex_substate_id, uint>& state) const {
    double res = 0;
    for (auto& filt: pfilters) {
        double rmult = 1.0;
        for (auto sus: filt.range()) {
            if (state[sus] < filt[sus].min or state[sus] > filt[sus].max) {
                rmult = 0;
                break;
            }
            uint nbLocked = filt[sus].min - preactants[sus];
            for (uint n = 0; n < preactants[sus]; ++n) {
                rmult *= static_cast<double>(state[sus] - n - nbLocked);
            }
            // divide by k!
            if (preactants[sus] > 1) {
                for (uint k = 2; k <= preactants[sus]; ++k) {
                    rmult /= static_cast<double>(k);
                }
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
    return ev2 != nullptr and ComplexLHSEventdef::isSame(ev) and preactants == ev2->reactants();
}

std::set<complex_substate_id> ComplexUpdateEventdef::getUpdSet() const {
    std::set<complex_substate_id> updset;
    for (auto upd: pupdates) {
        for (auto sus: upd.update.range()) {
            if (upd.update[sus] != 0) {
                updset.insert(sus);
            }
        }
    }
    return updset;
}

const util::strongid_vector<complex_substate_id, int>& ComplexUpdateEventdef::getUpdate(
    const util::strongid_vector<complex_substate_id, uint>& state,
    const rng::RNGptr& rng) const {
    uint ind = 0;
    if (pupdates.size() > 1) {
        ind = rng->get() % pupdates.size();
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
    // Should never reach there
    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

ComplexDeleteEventdef::ComplexDeleteEventdef(const model::ComplexDeleteEvent& ce,
                                             const Statedef& sd)
    : ComplexLHSEventdef::ComplexLHSEventdef(ce, sd) {}

double ComplexDeleteEventdef::rateMult(
    const util::strongid_vector<complex_substate_id, uint>& /*state*/) const {
    return 1.0;
}

std::set<complex_substate_id> ComplexDeleteEventdef::getUpdSet() const {
    std::set<complex_substate_id> updset;
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
    std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>> const& filts)
    const {
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

ComplexFilterChange::ComplexFilterChange(ComplexFilterChangeType tpe, complex_individual_id ind)
    : changeType(tpe)
    , stateInd(ind) {}

////////////////////////////////////////////////////////////////////////////////

ComplexFilter::ComplexFilter(
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& filts,
    complex_filter_id _id,
    const std::unordered_map<complex_individual_id, ComplexState>& states)
    : filterId(_id)
    , filters(filts)
    , pmatchAll(true) {
    AssertLog(filters.size() > 0);
    pDependsOnSus.container().resize(filters[0].size(), 0);
    for (auto& filt: filters) {
        for (auto sus: filt.range()) {
            pDependsOnSus[sus] |= filt[sus].min > 0 or
                                  filt[sus].max < model::COMPLEX_FILTER_MAX_VALUE;
            pmatchAll &= not pDependsOnSus[sus];
        }
    }

    for (auto stp: states) {
        toUpdate(stp.first);
    }
}

uint ComplexFilter::nbSus(
    const complex_substate_id& sus,
    const std::unordered_map<complex_individual_id, ComplexState>& states) const {
    uint cnt = 0;
    if (pmatchAll) {
        for (auto& stp: states) {
            cnt += stp.second[sus];
        }
    } else {
        for (auto ind: allMatches) {
            cnt += states.at(ind)[sus];
        }
    }
    return cnt;
}

void ComplexFilter::toUpdate(complex_individual_id ind) {
    ptoUpdate.insert(ind);
}

void ComplexFilter::reset() {
    ptoUpdate.clear();
    allMatches.clear();
    lastUpdates.clear();
}

void ComplexFilter::processUpdates(
    const std::unordered_map<complex_individual_id, ComplexState>& states) {
    if (pmatchAll) {
        for (auto ind: ptoUpdate) {
            lastUpdates.emplace_back(UPDMatch, ind);
        }
    } else {
        for (auto ind: ptoUpdate) {
            const auto st = states.find(ind);
            const bool match = (st != states.end()) ? computeMatch(st->second) : false;
            if (allMatches.find(ind) != allMatches.end()) {
                if (match) {
                    lastUpdates.emplace_back(UPDMatch, ind);
                } else {
                    allMatches.erase(ind);
                    lastUpdates.emplace_back(DELMatch, ind);
                }
            } else {
                if (match) {
                    allMatches.insert(ind);
                    lastUpdates.emplace_back(NEWMatch, ind);
                }
            }
        }
    }
    ptoUpdate.clear();
}

void ComplexFilter::clearLastUpdates() {
    lastUpdates.clear();
}

bool ComplexFilter::computeMatch(const ComplexState& state) const {
    if (pmatchAll) {
        return true;
    }

    for (auto& filt: filters) {
        bool match = true;
        for (auto sus: state.range()) {
            if (state[sus] < filt[sus].min or state[sus] > filt[sus].max) {
                match = false;
                break;
            }
        }
        if (match) {
            return true;
        }
    }
    return false;
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

ComplexState::ComplexState(const util::strongid_vector<complex_substate_id, uint>& state,
                           complex_individual_id _stateInd)
    : util::strongid_vector<complex_substate_id, uint>(state)
    , stateInd(_stateInd) {}

void ComplexState::checkpoint(std::ostream& cp_file) const {
    util::checkpoint(cp_file, stateInd);
    util::checkpoint(cp_file,
                     *dynamic_cast<const util::strongid_vector<complex_substate_id, uint>*>(this));
}

void ComplexState::restore(std::istream& cp_file) {
    util::restore(cp_file, stateInd);
    util::restore(cp_file, *dynamic_cast<util::strongid_vector<complex_substate_id, uint>*>(this));
}

}  // namespace steps::solver
