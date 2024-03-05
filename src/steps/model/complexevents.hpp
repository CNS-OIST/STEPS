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

#include "solver/fwd.hpp"
#include "util/common.hpp"

namespace steps::model {

enum ComplexLocation { COMP = 0, PATCH_IN = 1, PATCH_SURF = 2, PATCH_OUT = 3 };

extern const std::array<ComplexLocation, 3> AllPatchLocations;

inline constexpr uint COMPLEX_FILTER_MAX_VALUE = std::numeric_limits<uint>::max();

////////////////////////////////////////////////////////////////////////////////
/// Subunit state filter that holds the minimum and maximum number of subunits
/// that can be in a given state in order to match a filter.
struct SubunitStateFilter {
    SubunitStateFilter() = default;
    SubunitStateFilter(uint _min, uint _max)
        : min(_min)
        , max(_max) {}
    uint min{};
    uint max{};
};

inline bool operator==(const SubunitStateFilter& lhs, const SubunitStateFilter& rhs) {
    return lhs.min == rhs.min and lhs.max == rhs.max;
}

////////////////////////////////////////////////////////////////////////////////
/// A specific complex update that holds:
///     - the subunit states that are required for an update to be possible
///     - the update vector that represents the changes in subunit state numbers
struct ComplexUpdate {
    ComplexUpdate(std::vector<uint> req, std::vector<int> upd)
        : requirement(req)
        , update(upd) {}
    util::strongid_vector<solver::complex_substate_id, uint> requirement;
    util::strongid_vector<solver::complex_substate_id, int> update;
};

inline bool operator==(const ComplexUpdate& lhs, const ComplexUpdate& rhs) {
    return lhs.requirement == rhs.requirement and lhs.update == rhs.update;
}


////////////////////////////////////////////////////////////////////////////////
/// Base class for all complex events
/// A complex event describes what happens to a specific complex during a complex
/// reaction.
class ComplexEvent {
  public:
    ComplexEvent(std::string const& cmplxId);
    virtual ~ComplexEvent() = default;

    const std::string& complexId() const noexcept {
        return pcomplexId;
    };

  protected:
    const std::string pcomplexId;
};

////////////////////////////////////////////////////////////////////////////////
/// Complex event representing the creation of a new complex in a specific state
/// The corresponding complex only appears on the right hand side of a reaction.
class ComplexCreateEvent: public ComplexEvent {
  public:
    ComplexCreateEvent(std::string const& cmplxId, const std::vector<uint>& in);

    const std::vector<uint>& init() const noexcept {
        return pinit;
    }

  protected:
    const std::vector<uint> pinit;
};

////////////////////////////////////////////////////////////////////////////////
/// Base class for complex events that affect complexes present on the left hand
/// side of a reaction.
/// Contains a filter that determines which specific states match the reaction.
class ComplexLHSEvent: public ComplexEvent {
  public:
    ComplexLHSEvent(std::string const& cmplxId,
                    const std::vector<std::vector<SubunitStateFilter>>& filts);

    const std::vector<util::strongid_vector<solver::complex_substate_id, SubunitStateFilter>>&
    filters() const noexcept {
        return pfilters;
    }

  protected:
    std::vector<util::strongid_vector<solver::complex_substate_id, SubunitStateFilter>> pfilters;
};

////////////////////////////////////////////////////////////////////////////////
/// Complex event representing a modification in the state of a complex
/// It holds:
///     - a vector of reactants, that are used to indicate complex reactions that
///       directly involve subunit states
///     - a vector of possible updates
///     - a potential change of location
class ComplexUpdateEvent: public ComplexLHSEvent {
  public:
    ComplexUpdateEvent(std::string const& cmplxId,
                       const std::vector<std::vector<SubunitStateFilter>>& filts,
                       const std::vector<uint>& reac,
                       const std::vector<ComplexUpdate>& upd,
                       ComplexLocation destLoc);

    const std::vector<uint>& reactants() const noexcept {
        return preactants;
    }
    const std::vector<ComplexUpdate>& updates() const noexcept {
        return pupdates;
    }
    ComplexLocation destLoc() const noexcept {
        return pdestLoc;
    }

  protected:
    const std::vector<uint> preactants;
    const std::vector<ComplexUpdate> pupdates;
    const ComplexLocation pdestLoc;
};

////////////////////////////////////////////////////////////////////////////////
/// Complex event representing the deletion of a complex
class ComplexDeleteEvent: public ComplexLHSEvent {
  public:
    ComplexDeleteEvent(std::string const& cmplxId,
                       const std::vector<std::vector<SubunitStateFilter>>& filt);
};

}  // namespace steps::model

// END
