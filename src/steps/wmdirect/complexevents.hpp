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

// STL headers.
#include <array>
#include <unordered_map>
#include <unordered_set>

// STEPS headers.
#include "solver/complexeventsdef.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::wmdirect {

////////////////////////////////////////////////////////////////////////////////

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
class ComplexLHSCandidates {
  public:
    ComplexLHSCandidates(solver::complex_global_id cmplxId);

    template <typename T>
    void addEvent(const std::shared_ptr<const solver::ComplexLHSEventdef>& ev, T& cdef) {
        if (events[0] == nullptr) {
            events[0] = ev;
            filters[0] = cdef.GetFilter(complexIdx, ev->filters());
            nbEvents = 1;
        } else {
            AssertLog(events[1] == nullptr);
            events[1] = ev;
            filters[1] = cdef.GetFilter(complexIdx, ev->filters());
            nbEvents = 2;
            if (events[0]->sameReactants(events[1])) {
                sameReactants = true;
            }
        }
    }

    double rateMult(
        const std::unordered_map<solver::complex_individual_id, solver::ComplexState>& states);

    void _setRateMult(
        const std::unordered_map<solver::complex_individual_id, solver::ComplexState>& states);

    std::vector<
        std::pair<std::shared_ptr<const solver::ComplexLHSEventdef>, solver::complex_individual_id>>
    selectEvents(const rng::RNGptr& rng) const;

    void reset();

  protected:
    solver::complex_global_id complexIdx;
    uint nbEvents;
    bool sameReactants;
    std::array<std::shared_ptr<solver::ComplexFilter>, 2> filters;
    std::array<std::shared_ptr<const solver::ComplexLHSEventdef>, 2> events;
    std::array<std::unordered_map<solver::complex_individual_id, double>, 2> candidateMults;
    std::array<double, 2> totMult;
    std::array<double, 2> sumCommon;
    std::unordered_set<solver::complex_individual_id> commonCandidates;
    double rmult;
};

}  // namespace steps::wmdirect
