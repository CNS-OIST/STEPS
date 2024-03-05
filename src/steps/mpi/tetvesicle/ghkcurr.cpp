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

#include "mpi/tetvesicle/ghkcurr.hpp"

// Standard library & STL headers.
#include <cmath>
#include <iostream>
#include <vector>

// STEPS headers.
#include "math/constants.hpp"
#include "math/ghk.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

GHKcurr::GHKcurr(solver::GHKcurrdef* ghkdef, TriRDEF* tri)
    : pGHKcurrdef(ghkdef)
    , pTri(tri)
    , pEffFlux(true) {
    AssertLog(pGHKcurrdef != nullptr);
    AssertLog(pTri != nullptr);

    pType = KP_GHK;
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
    pEffFlux = true;  // TODO: come back to this and check if rate needs to be
                      // recalculated here
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::setupDeps() {
    AssertLog(pTri->getInHost());

    KProcPSet local;

    // The only concentration changes for a GHK current event are in the outer
    // and inner volume. The flux can involve movement of ion from either
    // compartment to the other- depnding on direction of flux

    TetRDEF* itet = pTri->iTet();
    TetRDEF* otet = pTri->oTet();
    AssertLog(itet != nullptr);

    if (itet->getHost() != pTri->getHost()) {
        std::ostringstream os;
        os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
        os << itet->idx() << " belong to different hosts.\n";
        NotImplErrLog(os.str());
    }

    // The global species of the ion
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();

    // First check KProcs in the inner tetrahedron
    uint nkprocs = itet->countKProcs();
    for (uint k = 0; k < nkprocs; k++) {
        if (itet->KProcDepSpecTet(k, itet, gidxion)) {
            local.insert(itet->getKProc(k));
        }
    }

    for (auto const& tri: itet->nexttris()) {
        if (tri == nullptr) {
            continue;
        }

        if (itet->getHost() != tri->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron ";
            os << itet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = tri->countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            if (tri->KProcDepSpecTet(sk, itet, gidxion)) {
                local.insert(tri->getKProc(sk));
            }
        }
    }

    if (otet != nullptr) {
        if (otet->getHost() != pTri->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron ";
            os << otet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        // Now check KProcs in the outer tetrahedron
        nkprocs = otet->countKProcs();
        for (uint k = 0; k < nkprocs; k++) {
            if (otet->KProcDepSpecTet(k, otet, gidxion)) {
                local.insert(otet->getKProc(k));
            }
        }

        for (auto const& tri: otet->nexttris()) {
            if (tri == nullptr) {
                continue;
            }

            if (otet->getHost() != tri->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << tri->idx() << " and its compartment tetrahedron ";
                os << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = tri->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                if (tri->KProcDepSpecTet(sk, otet, gidxion)) {
                    local.insert(tri->getKProc(sk));
                }
            }
        }
    }

    localUpdVec.assign(local.begin(), local.end());
}

////////////////////////////////////////////////////////////////////////////////

double GHKcurr::rate(TetVesicleRDEF* solver) {
    // The rate comes from the GHK flux equation. The flux equation
    // returns a single-channel current which must be converted to a rate,
    // remembering that flux can be positive or negative (bi-directional)
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();
    double voconc = pGHKcurrdef->voconc();

    // Get concentrations in Molar units: convert to Mol/m^3
    double iconc = (pTri->iTet()->conc(gidxion)) * 1.0e3;
    double oconc = 0.0;

    if (voconc < 0.0) {
        oconc = (pTri->oTet()->conc(gidxion)) * 1.0e3;
    } else {
        oconc = voconc * 1.0e3;
    }
    double v = solver->getTriV_(pTri->idx());
    double T = solver->getTemp();  // seems safe enough to call this simple API function here

    double flux = math::GHKcurrent(
        pGHKcurrdef->perm(), v + pGHKcurrdef->vshift(), pGHKcurrdef->valence(), T, iconc, oconc);

    // Note: For a positive flux, this could be an efflux of +ve cations,
    // or an influx of -ve anions. Need to check the valence.

    // Get the rate of ion flux, remembering valence may by other than 1.
    double rt = flux / (math::E_CHARGE * static_cast<double>(pGHKcurrdef->valence()));
    // Now a positive rate is always an efflux and a negative rate is an influx
    // Set positive or negative flux flag
    if (rt >= 0.0) {
        setEffFlux(true);
    } else {
        setEffFlux(false);
    }

    // Find the number of available channel states
    solver::Patchdef* pdef = pTri->patchdef();
    solver::ghkcurr_local_id ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());
    // Fetch the local index of the channelstate
    solver::spec_local_id cslidx = pdef->ghkcurr_chanstate(ghklidx);
    auto n = static_cast<double>(pTri->pools()[cslidx]);

    return fabs(rt) * n;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& GHKcurr::getRemoteUpdVec(int /*direction*/) const {
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& GHKcurr::getLocalUpdVec(int /*direction*/) const {
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::apply(const rng::RNGptr& /*rng*/,
                    double /*dt*/,
                    double /*simtime*/,
                    double period,
                    TetVesicleRDEF* /*solver*/) {
    TetRDEF* itet = pTri->iTet();
    TetRDEF* otet = pTri->oTet();

    solver::Compdef* innercdef = itet->compdef();
    solver::Compdef* outercdef = nullptr;
    if (otet != nullptr) {
        outercdef = otet->compdef();
    }

    // Fetch the global index of the ion and the valence
    const solver::spec_global_id gidxion = pGHKcurrdef->ion();
    const uint valence = pGHKcurrdef->valence();
    AssertLog(valence != 0);

    solver::Patchdef* pdef = pTri->patchdef();
    const solver::ghkcurr_local_id ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());

    const bool realflux = pGHKcurrdef->realflux();
    double voconc = pGHKcurrdef->voconc();

    solver::spec_local_id linneridx = innercdef->specG2L(gidxion);
    solver::spec_local_id louteridx;
    if (outercdef != nullptr) {
        louteridx = outercdef->specG2L(gidxion);
    }

    if (efflux()) {
        if (realflux) {
            if (itet->clamped(linneridx) == false) {
                itet->incCount(linneridx, -1, period);
            }
            if ((otet != nullptr) && voconc < 0.0) {
                if (otet->clamped(louteridx) == false) {
                    otet->incCount(louteridx, 1, period);
                }
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, valence);
    } else {
        if (realflux) {
            if (itet->clamped(linneridx) == false) {
                itet->incCount(linneridx, 1);
            }
            if ((otet != nullptr) && voconc < 0.0) {
                if (otet->clamped(louteridx) == false) {
                    otet->incCount(louteridx, -1);
                }
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, -valence);
    }

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

void GHKcurr::resetOccupancies() {
    pTri->resetPoolOccupancy();

    // Update inner tet pools.
    TetRDEF* itet = pTri->iTet();
    if (itet != nullptr) {
        itet->resetPoolOccupancy();
    }

    TetRDEF* otet = pTri->oTet();
    if (otet != nullptr) {
        otet->resetPoolOccupancy();
    }
}

}  // namespace steps::mpi::tetvesicle
