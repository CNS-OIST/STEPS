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

#pragma once

// STL headers.
#include <string>

// STEPS headers.
#include "fwd.hpp"

#include "model/complex.hpp"
#include "model/complexevents.hpp"
#include "rng/rng.hpp"
#include "util/checkpointing.hpp"
#include "util/error.hpp"
#include "util/vocabulary.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::dist {

enum ComplexEventType { UPDEvent = 0, DELEvent, CREEvent };

class ComplexEventdef {
  public:
    ComplexEventdef(const steps::model::ComplexEvent& ce, const Statedef& sd);
    virtual ~ComplexEventdef() = default;

    const model::complex_id& complexIdx() const noexcept {
        return pcomplexIdx;
    }

    virtual ComplexEventType type() const noexcept = 0;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const {
        return type() == ev->type() and pcomplexIdx == ev->complexIdx();
    }

    virtual std::set<model::complex_substate_id> getUpdSet() const = 0;

    virtual std::set<model::complex_substate_id> getDepSet() const = 0;

  protected:
    const model::complex_id pcomplexIdx;
};

class ComplexCreateEventdef: public ComplexEventdef {
  public:
    ComplexCreateEventdef(const steps::model::ComplexCreateEvent& ce, const Statedef& sd);

    const util::strongid_vector<model::complex_substate_id, uint>& init() const noexcept {
        return pinit;
    }

    virtual ComplexEventType type() const noexcept override {
        return CREEvent;
    }

    virtual std::set<model::complex_substate_id> getUpdSet() const override;

    virtual std::set<model::complex_substate_id> getDepSet() const override {
        return {};
    }

  protected:
    const util::strongid_vector<model::complex_substate_id, uint> pinit;
};

class ComplexLHSEventdef: public ComplexEventdef {
  public:
    ComplexLHSEventdef(const steps::model::ComplexLHSEvent& ce, const Statedef& sd);

    const std::vector<
        util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>&
    filters() const noexcept {
        return pfilters;
    }

    virtual double rateMult(
        const util::strongid_vector<model::complex_substate_id, uint>& state) const = 0;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual bool sameReactants(const std::shared_ptr<const ComplexEventdef>& /*ev*/) const {
        return false;
    }

    virtual bool hasNoReactants() const {
        return true;
    }

    /// Get the set of complex substates on which the event depends
    virtual std::set<model::complex_substate_id> getDepSet() const override;

  protected:
    std::vector<util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
        pfilters;
};

class ComplexUpdateEventdef: public ComplexLHSEventdef {
  public:
    ComplexUpdateEventdef(const steps::model::ComplexUpdateEvent& ce, const Statedef& sd);

    const util::strongid_vector<model::complex_substate_id, uint>& reactants() const noexcept {
        return preactants;
    }
    const std::vector<steps::model::ComplexUpdate>& updates() const noexcept {
        return pupdates;
    }

    double rateMult(
        const util::strongid_vector<model::complex_substate_id, uint>& state) const override;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual bool sameReactants(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual bool hasNoReactants() const override {
        return pNoReactants;
    }

    virtual ComplexEventType type() const noexcept override {
        return UPDEvent;
    }

    steps::model::Location destLoc() const noexcept {
        return pdestLoc;
    }

    virtual std::set<model::complex_substate_id> getUpdSet() const override;

    const util::strongid_vector<model::complex_substate_id, int>& getUpdate(
        const util::strongid_vector<model::complex_substate_id, uint>& state,
        rng::RNG& rng) const;

  protected:
    const util::strongid_vector<model::complex_substate_id, uint> preactants;
    const std::vector<steps::model::ComplexUpdate> pupdates;
    const steps::model::Location pdestLoc;
    bool pNoReactants;
};

class ComplexDeleteEventdef: public ComplexLHSEventdef {
  public:
    ComplexDeleteEventdef(const steps::model::ComplexDeleteEvent& ce, const Statedef& sd);

    double rateMult(
        const util::strongid_vector<model::complex_substate_id, uint>& state) const override;

    virtual ComplexEventType type() const noexcept override {
        return DELEvent;
    }

    virtual std::set<model::complex_substate_id> getUpdSet() const override;
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// State of a specific complex
/// Contains a pool of subunit states
template <typename Entity>
class ComplexState: public util::strongid_vector<model::complex_substate_id, uint> {
  public:
    ComplexState() = default;

    ComplexState(const util::strongid_vector<model::complex_substate_id, uint>& state,
                 model::complex_individual_id _stateInd,
                 Entity _loc)
        : util::strongid_vector<model::complex_substate_id, uint>(state)
        , stateInd(_stateInd)
        , loc(_loc) {}

    const model::complex_individual_id& ind() const noexcept {
        return stateInd;
    }

    const Entity& location() const noexcept {
        return loc;
    }

    void checkpoint(std::ostream& cp_file) const {
        util::checkpoint(cp_file, stateInd);
        util::checkpoint(
            cp_file,
            *dynamic_cast<const util::strongid_vector<model::complex_substate_id, uint>*>(this));
    }

    void restore(std::istream& cp_file) {
        util::restore(cp_file, stateInd);
        util::restore(cp_file,
                      *dynamic_cast<util::strongid_vector<model::complex_substate_id, uint>*>(
                          this));
    }

  protected:
    model::complex_individual_id stateInd;
    Entity loc;
};

////////////////////////////////////////////////////////////////////////////////

struct FilterHash {
    std::size_t operator()(
        std::vector<util::strongid_vector<model::complex_substate_id,
                                          steps::model::SubunitStateFilter>> const& filts) const;
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// Description of a complex filter
///
/// This class is used to describe a complex filter in the def layer while the
/// actual complex filters are not yet created.
struct ComplexFilterDescr {
    ComplexFilterDescr(const steps::model::ComplexFilterDescr& cfd, const Statedef& sd);
    model::complex_id complexId;
    std::vector<util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
        filters;
};

////////////////////////////////////////////////////////////////////////////////
/// Complex filter that keeps track of the specific complexes that match it
///
/// A filter is composed of several subfilters, the filter matches a complex state
/// if any of the subfilter matches the state. A subfilter contains a minimum and
/// maximum value for each subunit state; the subfilter matches a complex state if
/// all subunit states in the complex state are between the specified minimum and
/// maximum values.
///
/// Complex filters keep a set of specific complexes that match it. This set is
/// however not kept up to date at all times:
///     - Processes that modify complex states are responsible for signaling the
///       changes to ComplexFilter by calling the toUpdate() method.
///     - Processes that require the set of matches are responsible for calling
///       processUpdates(). They can then get only the changes in matches by
///       calling getLastUpdates(). Once these changes are processed by all
///       objects that need to process them, clearLastUpdates() should be called.
template <typename Entity>
class ComplexFilter {
  public:
    ComplexFilter(
        const std::vector<util::strongid_vector<model::complex_substate_id,
                                                steps::model::SubunitStateFilter>>& filts,
        model::complex_filter_id _id,
        const std::unordered_map<model::complex_individual_id, ComplexState<Entity>>& states)
        : filterId(_id)
        , filters(filts)
        , pmatchAll(true) {
        AssertLog(filters.size() > 0);
        pDependsOnSus.container().resize(filters[0].size(), 0);
        for (auto& filt: filters) {
            for (auto sus: filt.range()) {
                pDependsOnSus[sus] |= filt[sus].min > 0 or
                                      filt[sus].max < steps::model::COMPLEX_FILTER_MAX_VALUE;
                pmatchAll &= not pDependsOnSus[sus];
            }
        }

        for (auto stp: states) {
            toUpdate(stp.first);
        }
    }


    model::complex_filter_id id() const {
        return filterId;
    }

    std::optional<model::complex_filter_occupancy_id> occupancy_id() const {
        return occupancyId;
    }

    void setOccupancyId(model::complex_filter_occupancy_id id) {
        assert(not occupancyId);
        occupancyId = id;
    }

    bool dependsOnSus(model::complex_substate_id sus) const {
        return pDependsOnSus[sus] > 0;
    }

    inline bool matchAll() const {
        return pmatchAll;
    }

    uint nbMatches(const std::unordered_map<model::complex_individual_id, ComplexState<Entity>>&
                       states) const {
        return pmatchAll ? states.size() : allMatches.size();
    }

    inline void addMatch(model::complex_individual_id ind) {
        allMatches.insert(ind);
    }

    inline void removeMatch(model::complex_individual_id ind) {
        allMatches.erase(ind);
    }

    inline bool matches(model::complex_individual_id ind) const {
        return pmatchAll or allMatches.find(ind) != allMatches.end();
    }

    inline void toUpdate(model::complex_individual_id ind) {
        ptoUpdate.insert(ind);
    }

    void reset() {
        ptoUpdate.clear();
        allMatches.clear();
    }

    void processUpdates(
        const std::unordered_map<model::complex_individual_id, ComplexState<Entity>>& states) {
        if (not pmatchAll) {
            for (auto ind: ptoUpdate) {
                const auto st = states.find(ind);
                const bool match = (st != states.end()) ? computeMatch(st->second) : false;
                if (allMatches.find(ind) != allMatches.end()) {
                    if (not match) {
                        allMatches.erase(ind);
                    }
                } else {
                    if (match) {
                        allMatches.insert(ind);
                    }
                }
            }
        }
        ptoUpdate.clear();
    }

    bool isUpdated() const noexcept {
        return ptoUpdate.empty();
    }

    bool computeMatch(const ComplexState<Entity>& state) const {
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

  protected:
    const model::complex_filter_id filterId;
    const std::vector<
        util::strongid_vector<model::complex_substate_id, steps::model::SubunitStateFilter>>
        filters;

    std::optional<model::complex_filter_occupancy_id> occupancyId;

    std::set<model::complex_individual_id> ptoUpdate;
    std::unordered_set<model::complex_individual_id> allMatches;

    util::strongid_vector<model::complex_substate_id, uint> pDependsOnSus;
    bool pmatchAll;
};


////////////////////////////////////////////////////////////////////////////////
/// Candidate specific complexes of a given complex type
///
/// This class keeps track of the specific complexes that are candidates for
/// participating in a complex reaction. It can hold up to two complex filters
/// corresponding to two distinct complexes of the same type. Since some
/// specific complexes can match both filters but can only be involved in the
/// reaction one time, we need this class to keep track of the candidates that
/// are common to both filters and compute the rate multiplier accordingly.
/// In addition, this class is also responsible for sampling the specific
/// complexes that will be used during a complex reaction.
template <typename Entity>
class ComplexLHSCandidates {
  public:
    ComplexLHSCandidates(model::complex_id cmplxId, Entity ent);

    template <typename T>
    void addEvent(const std::shared_ptr<const ComplexLHSEventdef>& ev, const T& mol_ent) {
        if (events[0] == nullptr) {
            events[0] = ev;
            filters[0] = mol_ent.GetComplexFilter(complexIdx_, ev->filters());
            nbEvents = 1;
        } else {
            AssertLog(events[1] == nullptr);
            events[1] = ev;
            filters[1] = mol_ent.GetComplexFilter(complexIdx_, ev->filters());
            nbEvents = 2;
            sameReactants = events[0]->sameReactants(events[1]);
        }
    }

    double rateMult(const EntityMolecules<Entity>& entmols) const {
        const auto& states = entmols.complexStates(complexIdx_);
        for (uint i = 0; i < nbEvents; ++i) {
            filters[i]->processUpdates(states);
        }
        double diag = 0;
        double totMult0 = 0, totMult1 = 0;
        double sumCommon0 = 0, sumCommon1 = 0;

        auto [bgn, end] = entmols.complexLocations(complexIdx_).equal_range(entity);
        if (events[0]->hasNoReactants() and (nbEvents == 1 or events[1]->hasNoReactants())) {
            // No need to compute rate multipliers per match if both events have no reactants
            for (auto it = bgn; it != end; ++it) {
                if (filters[0]->matches(it->second)) {
                    totMult0 += 1;
                    if (nbEvents > 1 and filters[1]->matches(it->second)) {
                        totMult1 += 1;
                        sumCommon0 += 1;
                        sumCommon1 += 1;
                        diag += 1;
                    }
                } else if (nbEvents > 1 and filters[1]->matches(it->second)) {
                    totMult1 += 1;
                }
            }
        } else {
            for (auto it = bgn; it != end; ++it) {
                if (filters[0]->matches(it->second)) {
                    auto stit = states.find(it->second);
                    if (stit != states.end()) {
                        double rm0 = events[0]->rateMult(stit->second);
                        totMult0 += rm0;
                        // If other filter matches
                        if (nbEvents > 1 and filters[1]->matches(it->second)) {
                            double rm1 = events[1]->rateMult(stit->second);
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
                        totMult1 += rm1;
                    }
                }
            }
        }

        if (nbEvents > 1) {
            if (sameReactants) {
                return totMult0 * totMult1 - (sumCommon0 * sumCommon1 + diag) / 2.0;
            } else {
                return totMult0 * totMult1 - diag;
            }
        } else {
            return totMult0;
        }
    }

    std::vector<std::pair<std::shared_ptr<const ComplexLHSEventdef>, model::complex_individual_id>>
    selectEvents(const EntityMolecules<Entity>& entmols, rng::RNG& rng) const;

    model::complex_id complexIdx() const noexcept {
        return complexIdx_;
    }

  protected:
    model::complex_id complexIdx_;
    uint nbEvents;
    bool sameReactants;
    std::array<std::shared_ptr<ComplexFilter<Entity>>, 2> filters;
    std::array<std::shared_ptr<const ComplexLHSEventdef>, 2> events;
    Entity entity;
};

}  // namespace steps::dist
