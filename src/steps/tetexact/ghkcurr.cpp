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

// STEPS headers.
#include "ghkcurr.hpp"
#include "math/ghk.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "tri.hpp"
#include "wmvol.hpp"

#include "util/error.hpp"
#include <easylogging++.h>

#include "util/checkpointing.hpp"

namespace steps::tetexact {

GHKcurr::GHKcurr(solver::GHKcurrdef* ghkdef, Tri* tri)
    : pGHKcurrdef(ghkdef)
    , pTri(tri)
    , pEffFlux(true) {
    AssertLog(pGHKcurrdef != nullptr);
    AssertLog(pTri != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr::~GHKcurr() = default;

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pEffFlux);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::restore(std::fstream& cp_file) {
    util::restore(cp_file, pEffFlux);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::reset() {
    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    setActive(true);
    pEffFlux = true;  // TODO: come back to this and check if rate needs to be recalculated here
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setupDeps() {
    KProcPSet updset;

    // The only concentration changes for a GHK current event are in the outer
    // and inner volume. The flux can involve movement of ion from either
    // compartment to the other- depnding on direction of flux

    WmVol* itet = pTri->iTet();
    WmVol* otet = pTri->oTet();
    AssertLog(itet != nullptr);
    // The global species of the ion
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();

    // First check KProcs in the inner tetrahedron
    for (auto const& k: itet->kprocs()) {
        if (k->depSpecTet(gidxion, itet)) {
            updset.insert(k);
        }
    }

    for (auto const& tri: itet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }

        for (auto const& k: tri->kprocs()) {
            if (k->depSpecTet(gidxion, itet)) {
                updset.insert(k);
            }
        }
    }

    if (otet != nullptr) {
        // Now check KProcs in the outer tetrahedron
        for (auto const& k: otet->kprocs()) {
            if (k->depSpecTet(gidxion, otet)) {
                updset.insert(k);
            }
        }

        for (auto const& tri: otet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }

            for (auto const& k: tri->kprocs()) {
                if (k->depSpecTet(gidxion, otet)) {
                    updset.insert(k);
                }
            }
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

bool GHKcurr::depSpecTet(solver::spec_global_id gidx, WmVol* tet) {
    // Two things here. Any changes of the count of the ion in the inner or
    // outer tetrahedron affect the single-channel rate. This will be dealt with
    // together for now in rate()
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp
    if (tet == pTri->iTet()) {
        return (pGHKcurrdef->dep_v(gidx) != solver::DEP_NONE);
    } else if (tet == pTri->oTet()) {
        if (pGHKcurrdef->voconc() < 0.0) {
            return (pGHKcurrdef->dep_v(gidx) != solver::DEP_NONE);
        } else {
            return false;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool GHKcurr::depSpecTri(solver::spec_global_id gidx, Tri* triangle) {
    if (triangle != pTri) {
        return false;
    }
    return (pGHKcurrdef->dep(gidx) != solver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::rate(steps::tetexact::Tetexact* solver) {
    // The rate comes from the GHK flux equation. The flux equation
    // returns a single-channel current which must be converted to a rate,
    // remembering that flux can be positive or negative (bi-directional)
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();
    double voconc = pGHKcurrdef->voconc();
    double oconc = 0.0;

    // Get concentrations in Molar units: convert to Mol/m^3
    double iconc = (pTri->iTet()->conc(gidxion)) * 1.0e3;

    if (voconc < 0.0) {
        oconc = (pTri->oTet()->conc(gidxion)) * 1.0e3;
    } else {
        oconc = voconc * 1.0e3;
    }

    double v = solver->getTriV(pTri->idx());
    double T = solver->getTemp();

    double flux = math::GHKcurrent(
        pGHKcurrdef->perm(), v + pGHKcurrdef->vshift(), pGHKcurrdef->valence(), T, iconc, oconc);

    // Note: For a positive flux, this could be an efflux of +ve cations,
    // or an influx of -ve anions. Need to check the valence.

    // Get the rate of ion flux, remembering valence may by other than 1.
    double rt = flux / (math::E_CHARGE * static_cast<double>(pGHKcurrdef->valence()));
    // Now a positive rate is always an efflux and a negative rate is an influx
    // Set positive or negative flux flag
    setEffFlux(rt >= 0.0);

    // Find the number of available channel states
    solver::Patchdef* pdef = pTri->patchdef();
    solver::ghkcurr_local_id ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());
    // Fetch the local index of the channelstate
    solver::spec_local_id cslidx = pdef->ghkcurr_chanstate(ghklidx);
    auto n = static_cast<double>(pTri->pools()[cslidx]);

    return fabs(rt) * n;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& GHKcurr::apply(const rng::RNGptr& /*rng*/,
                                          double /*dt*/,
                                          double /*simtime*/) {
    WmVol* itet = pTri->iTet();
    WmVol* otet = pTri->oTet();

    solver::Compdef* innercdef = itet->compdef();
    solver::Compdef* outercdef = nullptr;
    if (otet != nullptr) {
        outercdef = otet->compdef();
    }

    // Fetch the global index of the ion and the valence
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();
    const auto valence = pGHKcurrdef->valence();

    solver::Patchdef* pdef = pTri->patchdef();
    const solver::ghkcurr_local_id ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());

    AssertLog(valence != 0);  // TODO: get rid of this when tested

    const bool realflux = pGHKcurrdef->realflux();
    double voconc = pGHKcurrdef->voconc();

    solver::spec_local_id linneridx = innercdef->specG2L(gidxion);

    solver::spec_local_id louteridx;
    if (outercdef != nullptr) {
        louteridx = outercdef->specG2L(gidxion);
    }

    if (efflux()) {
        if (realflux) {
            if (!itet->clamped(linneridx)) {
                itet->incCount(linneridx, -1);
            }
            if ((otet != nullptr) && voconc < 0.0) {
                if (!otet->clamped(louteridx)) {
                    otet->incCount(louteridx, 1);
                }
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, valence);
    } else {
        if (realflux) {
            if (!itet->clamped(linneridx)) {
                itet->incCount(linneridx, 1);
            }
            if ((otet != nullptr) && voconc < 0.0) {
                if (!otet->clamped(louteridx)) {
                    otet->incCount(louteridx, -1);
                }
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, -valence);
    }

    rExtent++;

    return pUpdVec;
}

}  // namespace steps::tetexact
