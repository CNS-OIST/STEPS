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

#pragma once

#include <set>

#include "model/complex.hpp"
#include "model/complexevents.hpp"
#include "solver/statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

enum ComplexEventType { UPDEvent = 0, DELEvent, CREEvent };

class ComplexEventdef {
  public:
    ComplexEventdef(const model::ComplexEvent& ce, const Statedef& sd);
    virtual ~ComplexEventdef() = default;

    const complex_global_id& complexIdx() const noexcept {
        return pcomplexIdx;
    }

    virtual ComplexEventType type() const noexcept = 0;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const {
        return type() == ev->type() and pcomplexIdx == ev->complexIdx();
    }

    virtual std::set<complex_substate_id> getUpdSet() const = 0;

    virtual std::set<complex_substate_id> getDepSet() const = 0;

  protected:
    const complex_global_id pcomplexIdx;
};

class ComplexCreateEventdef: public ComplexEventdef {
  public:
    ComplexCreateEventdef(const model::ComplexCreateEvent& ce, const Statedef& sd);

    const util::strongid_vector<complex_substate_id, uint>& init() const noexcept {
        return pinit;
    }

    virtual ComplexEventType type() const noexcept override {
        return CREEvent;
    }

    virtual std::set<complex_substate_id> getUpdSet() const override;

    virtual std::set<complex_substate_id> getDepSet() const override {
        return {};
    }

  protected:
    const util::strongid_vector<complex_substate_id, uint> pinit;
};

class ComplexLHSEventdef: public ComplexEventdef {
  public:
    ComplexLHSEventdef(const model::ComplexLHSEvent& ce, const Statedef& sd);

    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>&
    filters() const noexcept {
        return pfilters;
    }

    virtual double rateMult(
        const util::strongid_vector<complex_substate_id, uint>& state) const = 0;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual bool sameReactants(const std::shared_ptr<const ComplexEventdef>& /*ev*/) const {
        return false;
    }

    virtual std::set<complex_substate_id> getDepSet() const override;

  protected:
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>&
        pfilters;
};

class ComplexUpdateEventdef: public ComplexLHSEventdef {
  public:
    ComplexUpdateEventdef(const model::ComplexUpdateEvent& ce, const Statedef& sd);

    const util::strongid_vector<complex_substate_id, uint>& reactants() const noexcept {
        return preactants;
    }
    const std::vector<model::ComplexUpdate>& updates() const noexcept {
        return pupdates;
    }

    double rateMult(const util::strongid_vector<complex_substate_id, uint>& state) const override;

    virtual bool isSame(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual bool sameReactants(const std::shared_ptr<const ComplexEventdef>& ev) const override;

    virtual ComplexEventType type() const noexcept override {
        return UPDEvent;
    }

    model::ComplexLocation destLoc() const noexcept {
        return pdestLoc;
    }

    virtual std::set<complex_substate_id> getUpdSet() const override;

    const util::strongid_vector<complex_substate_id, int>& getUpdate(
        const util::strongid_vector<complex_substate_id, uint>& state,
        const rng::RNGptr& rng) const;

  protected:
    const util::strongid_vector<complex_substate_id, uint> preactants;
    const std::vector<model::ComplexUpdate> pupdates;
    const model::ComplexLocation pdestLoc;
};

class ComplexDeleteEventdef: public ComplexLHSEventdef {
  public:
    ComplexDeleteEventdef(const model::ComplexDeleteEvent& ce, const Statedef& sd);

    double rateMult(const util::strongid_vector<complex_substate_id, uint>& state) const override;

    virtual ComplexEventType type() const noexcept override {
        return DELEvent;
    }

    virtual std::set<complex_substate_id> getUpdSet() const override;
};

////////////////////////////////////////////////////////////////////////////////

////////////////////////////////////////////////////////////////////////////////
/// State of a specific complex
/// Contains a pool of subunit states
class ComplexState: public util::strongid_vector<complex_substate_id, uint> {
  public:
    ComplexState() = default;

    ComplexState(const util::strongid_vector<complex_substate_id, uint>& state,
                 complex_individual_id _stateInd);

    const complex_individual_id& ind() const {
        return stateInd;
    }

    void checkpoint(std::ostream& cp_file) const;

    void restore(std::istream& cp_file);

  protected:
    complex_individual_id stateInd;
};

////////////////////////////////////////////////////////////////////////////////

struct FilterHash {
    std::size_t operator()(
        std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>> const&
            filts) const;
};

////////////////////////////////////////////////////////////////////////////////

enum ComplexFilterChangeType { UPDMatch = 0, DELMatch, NEWMatch };

////////////////////////////////////////////////////////////////////////////////
/// Utility struct for specifying a type of change in a match and the
/// corresponding specific complex.
struct ComplexFilterChange {
    ComplexFilterChange(ComplexFilterChangeType tpe, complex_individual_id ind);

    ComplexFilterChangeType changeType;
    complex_individual_id stateInd;
};

////////////////////////////////////////////////////////////////////////////////

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
class ComplexFilter {
  public:
    enum MatchType { MATCHFalse = 0, MATCHTrue, MATCHUnknown };

    ComplexFilter(const std::vector<
                      util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>& filts,
                  complex_filter_id _id,
                  const std::unordered_map<complex_individual_id, ComplexState>& states);

    complex_filter_id id() const {
        return filterId;
    }

    bool dependsOnSus(complex_substate_id sus) const {
        return pDependsOnSus[sus];
    }

    bool matchAll() const {
        return pmatchAll;
    }

    uint nbMatches(const std::unordered_map<complex_individual_id, ComplexState>& states) const {
        return pmatchAll ? states.size() : allMatches.size();
    }

    const std::vector<ComplexFilterChange>& getLastUpdates() const {
        return lastUpdates;
    }

    uint nbSus(const complex_substate_id& sus,
               const std::unordered_map<complex_individual_id, ComplexState>& states) const;

    bool matches(complex_individual_id ind) const {
        return allMatches.find(ind) != allMatches.end();
    }

    void toUpdate(complex_individual_id);

    void reset();

    void processUpdates(const std::unordered_map<complex_individual_id, ComplexState>& states);

    void clearLastUpdates();

  protected:
    bool computeMatch(const ComplexState& state) const;

    const complex_filter_id filterId;
    const std::vector<util::strongid_vector<complex_substate_id, model::SubunitStateFilter>>
        filters;

    std::set<complex_individual_id> ptoUpdate;
    std::unordered_set<complex_individual_id> allMatches;

    std::vector<ComplexFilterChange> lastUpdates;

    util::strongid_vector<complex_substate_id, uint> pDependsOnSus;
    bool pmatchAll;
};

}  // namespace steps::solver
