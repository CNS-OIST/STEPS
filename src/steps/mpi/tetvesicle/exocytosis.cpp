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

#include "mpi/tetvesicle/exocytosis.hpp"

// Standard library & STL headers.
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

Exocytosis::Exocytosis(solver::Exocytosisdef* exodef, TetRDEF* tet)
    : pExocytosisdef(exodef)
    , pTet(tet)
    , pRate_zero(true) {
    AssertLog(pExocytosisdef != nullptr);
    AssertLog(pTet != nullptr);

    pType = KP_EXO;

    // If this is added to all tets then it's possible the rate could always
    // be zero. Dumb, but possible. Look at removing this later.
    for (auto const& tetnexttri: pTet->nexttris()) {
        if (tetnexttri == nullptr) {
            continue;
        }
        // We got membrane
        pRate_zero = false;
        break;
    }
}

////////////////////////////////////////////////////////////////////////////////

Exocytosis::~Exocytosis() = default;

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::checkpoint(std::fstream& cp_file) {
    // NOTE:
    //   pRate_zero, localUpdVec, remoteUpdVec are part of setup process
    //   Current plan is to regenrate pRate_per_ves and pTotal_rate after data restore
    //   as the 3rd part of the restoration. This should include a check that
    //  crData.rate == pTotal_rate (or just not necessary to checkpoint crData at all?)
    util::checkpoint(cp_file, pRate_zero);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::restore(std::fstream& cp_file) {
    util::compare(cp_file, pRate_zero);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    resetExtent();
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::setupDeps() {
    // Mark all exocytosis events for a potential update (including this one) since they may also
    // reference this same vesproxy, which will efectively be frozen upon execution of this event
    for (auto const& kproc: pTet->kprocs()) {
        if (kproc->getType() == KP_EXO) {
            localUpdVec.emplace_back(kproc);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

double Exocytosis::rate(TetVesicleRDEF* /*solver*/) {
    if (pRate_zero || inactive()) {
        return 0.0;
    }

    // This is an all-or-nothing type reaction. Need to check whether the
    // dependency species are present where they are supposed to be,
    // on the surface

    // Exodef object stores the 'LHS' information by global species index.
    // Let's just use that instead of the normal local STEPS indices
    // through the parent def object

    // This has two criteria now - the species have to be available and also there
    // must be some surface patch.

    auto const& vesproxyrefs = pTet->getVesProxyrefs();

    pRate_per_ves.clear();
    pTotal_rate = 0.0;

    // We don't bother with local indices here- just use the global ones
    const auto& lhs_v_vec = exodef()->lhs_V();

    for (auto const& vesproxy_mapit: vesproxyrefs) {
        auto& vesproxy_uid = vesproxy_mapit.first;
        auto& vesproxy = vesproxy_mapit.second;

        AssertLog(pRate_per_ves.find(vesproxy_uid) == pRate_per_ves.end());

        if (vesproxy->hasLinkSpec() || vesproxy->exoApplied().valid()) {
            // Zero rate if vesproxy has linkspec or exocytosis has already been applied
            pRate_per_ves[vesproxy_uid] = 0.0;
            continue;
        }

        if (exodef()->vdepssize() == 0) {
            // we know we have membrane (if pRate_zero==false) so this is simple
            pRate_per_ves[vesproxy_uid] = kcst();
            pTotal_rate += kcst();
            continue;
        } else  // not strictly necessary because we have the continue, but for
                // safety (in case that's ever changed)
        {
            for (auto sg: lhs_v_vec.range()) {
                uint lhs = lhs_v_vec[sg];
                if (lhs == 0) {
                    continue;
                }
                uint tetspeccnt = vesproxy->getSpecCount_V(sg);

                if (lhs > tetspeccnt) {
                    // The required species are not available
                    pRate_per_ves[vesproxy_uid] = 0.0;
                } else {
                    // The required species ARE available
                    pRate_per_ves[vesproxy_uid] = kcst();
                    pTotal_rate += kcst();
                }
            }
        }
    }

    return pTotal_rate;
}

////////////////////////////////////////////////////////////////////////////////

void Exocytosis::apply(const rng::RNGptr& rng,
                       double /*dt*/,
                       double /*simtime*/,
                       double /*period*/,
                       TetVesicleRDEF* /*solver*/) {
    AssertLog(pTotal_rate != 0);

    double selector = rng->getUnfII() * pTotal_rate;

    // Select by rate
    solver::vesicle_individual_id vesproxy_uid;
    double accum = 0.0;
    for (auto const& vesr: pRate_per_ves) {
        vesproxy_uid = vesr.first;
        accum += vesr.second;
        if (selector < accum) {
            break;
        }
    }

    // Each exocytosis event should only be selected once
    AssertLog(vesproxy_uid.valid());

    VesProxy* vesproxy = pTet->getVesProxyref(vesproxy_uid);

    if (vesproxy->exoApplied().valid()) {
        ProgErrLog("Vesicle has already been marked for exocytosis");
    }

    // Now set a flag to exocytose this vesicle. Freeze it in its current state-
    // it can undergo no further interactions then if marked for exocytosis (not
    // strictly necessary but just a bit more realistic, consistent)
    vesproxy->applyExo(exodef()->gidx());

    // That's it.
}

}  // namespace steps::mpi::tetvesicle
