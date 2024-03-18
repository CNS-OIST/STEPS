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

#include "wmdirect/complexevents.hpp"
#include "rng/rng.hpp"
#include "solver/compdef.hpp"
#include "util/error.hpp"

namespace steps::wmdirect {

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////

ComplexLHSCandidates::ComplexLHSCandidates(solver::complex_global_id cmplxId)
    : complexIdx(cmplxId)
    , nbEvents(0)
    , sameReactants(false)
    , filters{nullptr, nullptr}
    , events{nullptr, nullptr}
    , totMult({0.0, 0.0})
    , sumCommon({0.0, 0.0})
    , rmult(0) {}

double ComplexLHSCandidates::rateMult(
    const std::unordered_map<solver::complex_individual_id, solver::ComplexState>& states) {
    for (uint i = 0; i < nbEvents; ++i) {
        filters[i]->processUpdates(states);
    }
    for (uint i = 0; i < nbEvents; ++i) {
        for (const auto& ev: filters[i]->getLastUpdates()) {
            switch (ev.changeType) {
            case solver::UPDMatch: {
                const auto& state = states.at(ev.stateInd);
                const double newRate = events[i]->rateMult(state);
                const auto it = candidateMults[i].find(ev.stateInd);
                AssertLog(it != candidateMults[i].end());
                // Update the rate
                totMult[i] += newRate - it->second;
                it->second = newRate;
                break;
            }
            case solver::NEWMatch: {
                const auto& state = states.at(ev.stateInd);
                const double newRate = events[i]->rateMult(state);
                candidateMults[i].emplace(ev.stateInd, newRate);
                totMult[i] += newRate;

                if (nbEvents > 1 and filters[(i + 1) % nbEvents]->matches(ev.stateInd)) {
                    commonCandidates.insert(ev.stateInd);
                }
                break;
            }
            case solver::DELMatch: {
                const auto it = candidateMults[i].find(ev.stateInd);
                if (it != candidateMults[i].end()) {
                    totMult[i] -= it->second;
                    candidateMults[i].erase(it);
                }
                const auto it2 = commonCandidates.find(ev.stateInd);
                if (it2 != commonCandidates.end()) {
                    commonCandidates.erase(it2);
                }
                break;
            }
            default:
                AssertLog(false);
            }
        }
    }

    // Update ratemult
    if (nbEvents == 1) {
        rmult = totMult[0];
    } else {
        _setRateMult(states);
    }

    return rmult;
}

void ComplexLHSCandidates::_setRateMult(
    const std::unordered_map<solver::complex_individual_id, solver::ComplexState>& /*states*/) {
    // Correct the double counting, if needed
    double diag = 0.0;
    sumCommon[0] = 0.0;
    sumCommon[1] = 0.0;
    for (const auto& cand: commonCandidates) {
        const double rm1 = candidateMults[0][cand];
        const double rm2 = candidateMults[1][cand];
        diag += rm1 * rm2;
        sumCommon[0] += rm1;
        sumCommon[1] += rm2;
    }
    if (sameReactants) {
        rmult = totMult[0] * totMult[1] - (sumCommon[0] * sumCommon[1] + diag) / 2.0;
    } else {
        rmult = totMult[0] * totMult[1] - diag;
    }
}

std::vector<
    std::pair<std::shared_ptr<const solver::ComplexLHSEventdef>, solver::complex_individual_id>>
ComplexLHSCandidates::selectEvents(const steps::rng::RNGptr& rng) const {
    std::vector<
        std::pair<std::shared_ptr<const solver::ComplexLHSEventdef>, solver::complex_individual_id>>
        selected;

    // Start with the event that has the least choice
    std::array<uint, 2> ord{};
    if (nbEvents > 1 and candidateMults[0].size() > candidateMults[1].size()) {
        ord = {1, 0};
    } else {
        ord = {0, 1};
    }

    if (nbEvents > 1) {
        uint i = ord[0];
        uint j = ord[1];
        if (sameReactants) {
            // If the reactants are identical
            double p = rng->getUnfIE() * rmult;
            double tmp = 0.0;
            bool inCommon = false;
            // Select first complex state
            for (const auto& cand: candidateMults[i]) {
                inCommon = commonCandidates.find(cand.first) != commonCandidates.end();
                if (inCommon) {
                    tmp += cand.second *
                           (totMult[j] - (sumCommon[j] + candidateMults[j].at(cand.first)) / 2.0);
                } else {
                    tmp += cand.second * totMult[j];
                }
                if (tmp >= p) {
                    selected.emplace_back(events[i], cand.first);
                    break;
                }
            }
            // Select second complex state
            if (inCommon) {
                solver::complex_individual_id ind = selected.back().second;
                p = rng->getUnfIE() *
                    (totMult[j] - (sumCommon[j] + candidateMults[j].at(ind)) / 2.0);
                tmp = 0.0;
                for (const auto& cand: candidateMults[j]) {
                    if (cand.first == ind) {
                        continue;
                    }
                    if (commonCandidates.find(cand.first) != commonCandidates.end()) {
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
                p = rng->getUnfIE() * totMult[j];
                tmp = 0.0;
                for (const auto& cand: candidateMults[j]) {
                    tmp += cand.second;
                    if (tmp >= p) {
                        selected.emplace_back(events[j], cand.first);
                        break;
                    }
                }
            }
        } else {
            // If the reactants are different
            double p = rng->getUnfIE() * rmult;
            double tmp = 0.0;
            bool inCommon = false;
            // Select first complex state
            for (const auto& cand: candidateMults[i]) {
                inCommon = commonCandidates.find(cand.first) != commonCandidates.end();
                if (inCommon) {
                    tmp += cand.second * (totMult[j] - candidateMults[j].at(cand.first));
                } else {
                    tmp += cand.second * totMult[j];
                }
                if (tmp >= p) {
                    selected.emplace_back(events[i], cand.first);
                    break;
                }
            }
            // Select second complex state
            if (inCommon) {
                solver::complex_individual_id ind = selected.back().second;
                p = rng->getUnfIE() * (totMult[j] - candidateMults[j].at(ind));
                tmp = 0.0;
                for (const auto& cand: candidateMults[j]) {
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
                p = rng->getUnfIE() * totMult[j];
                tmp = 0.0;
                for (const auto& cand: candidateMults[j]) {
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
        const double p = rng->getUnfIE() * rmult;
        double tmp = 0.0;
        for (const auto& cand: candidateMults[0]) {
            tmp += cand.second;
            if (tmp >= p) {
                selected.emplace_back(events[0], cand.first);
                break;
            }
        }
    }

    return selected;
}

void ComplexLHSCandidates::reset() {
    commonCandidates.clear();
    for (uint i = 0; i < nbEvents; ++i) {
        totMult[i] = 0.0;
        filters[i]->reset();
        candidateMults[i].clear();
    }
}

}  // namespace steps::wmdirect
