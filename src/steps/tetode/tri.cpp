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
#include "tri.hpp"

#include "math/constants.hpp"
#include "math/ghk.hpp"
#include "solver/ghkcurrdef.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "tet.hpp"
#include "tetode.hpp"

#include "util/checkpointing.hpp"

namespace steps::tetode {

Tri::Tri(triangle_global_id idx,
         solver::Patchdef* patchdef,
         double area,
         double l0,
         double l1,
         double l2,
         double d0,
         double d1,
         double d2,
         tetrahedron_global_id tetinner,
         tetrahedron_global_id tetouter,
         triangle_global_id tri0,
         triangle_global_id tri1,
         triangle_global_id tri2)
    : pIdx(idx)
    , pPatchdef(patchdef)
    , pNextTri()
    , pArea(area)
    , pLengths()
    , pDist()

{
    AssertLog(pPatchdef != nullptr);
    AssertLog(pArea > 0.0);

    AssertLog(l0 > 0.0 && l1 > 0.0 && l2 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0);

    pTets[0] = tetinner;
    pTets[1] = tetouter;

    pTris[0] = tri0;
    pTris[1] = tri1;
    pTris[2] = tri2;

    pNextTri[0] = nullptr;
    pNextTri[1] = nullptr;
    pNextTri[2] = nullptr;

    pLengths[0] = l0;
    pLengths[1] = l1;
    pLengths[2] = l2;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
}

////////////////////////////////////////////////////////////////////////////////

Tri::~Tri() = default;

////////////////////////////////////////////////////////////////////////////////

void Tri::setInnerTet(Tet* t) {
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setOuterTet(Tet* t) {
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setNextTri(uint i, Tri* t) {
    AssertLog(i <= 2);

    pNextTri[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pERev);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::restore(std::fstream& cp_file) {
    util::restore(cp_file, pERev);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setOCerev(solver::ohmiccurr_local_id oclidx, double erev) {
    AssertLog(oclidx < patchdef()->countOhmicCurrs());
    pERev[oclidx] = erev;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getOCerev(solver::ohmiccurr_local_id oclidx) const {
    auto it = pERev.find(oclidx);
    if (it != pERev.end()) {
        return it->second;
    } else {
        return patchdef()->ohmiccurrdef(oclidx)->getERev();
    }
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getOhmicI(double v, steps::tetode::TetODE* solver) const {
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (auto i: solver::ohmiccurr_local_id::range(nocs)) {
        solver::OhmicCurrdef* ocdef = patchdef()->ohmiccurrdef(i);

        // Now need to get the states from TetODE object, and remember to convert local indices to
        // the global ones it needs
        solver::spec_global_id spec_gidx = patchdef()->specL2G(patchdef()->ohmiccurr_chanstate(i));

        double n = solver->_getTriSpecCount(pIdx, spec_gidx);

        current += (n * ocdef->getG()) * (v - getOCerev(i));
    }

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getGHKI(double v, double dt, steps::tetode::TetODE* solver) const {
    double current = 0.0;

    uint nghkcurrs = patchdef()->countGHKcurrs();
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        solver::GHKcurrdef* ghkdef = patchdef()->ghkcurrdef(i);

        // The rate comes from the GHK flux equation. The flux equation
        // returns a single-channel current which must be converted to a rate,
        // remembering that flux can be positive or negative (bi-directional)
        const solver::spec_global_id gidxion = ghkdef->ion();
        double voconc = ghkdef->voconc();

        // Get concentrations in Molar units: convert to Mol/m^3
        double iconc = solver->_getTetSpecConc(iTet()->idx(), gidxion) * 1.0e3;
        double oconc = 0.0;

        if (voconc < 0.0) {
            oconc = solver->_getTetSpecConc(oTet()->idx(), gidxion) * 1.0e3;
        } else {
            oconc = voconc * 1.0e3;
        }

        double T = solver->getTemp();

        double flux_per_channel = math::GHKcurrent(
            ghkdef->perm(), v + ghkdef->vshift(), ghkdef->valence(), T, iconc, oconc);

        // Fetch global index of channel state
        solver::spec_global_id cs_gidx = ghkdef->chanstate();

        double flux = solver->_getTriSpecCount(idx(), cs_gidx) * flux_per_channel;

        current += flux;

        // Now got to move ions
        if (ghkdef->realflux()) {
            // Note: For a positive flux, this could be an efflux of +ve cations,
            // or an influx of -ve anions. Need to check the valence.

            // Get the rate of ion flux, remembering valence may by other than 1.
            double rt = flux / (math::E_CHARGE * static_cast<double>(ghkdef->valence()));
            // Now a positive rate is always an efflux and a negative rate is an influx

            // rt is number of ions per second; positive is an efflux and a negative is an influx
            double count = rt * dt;

            if (voconc < 0.0) {
                solver->_setTetSpecCount(oTet()->idx(),
                                         gidxion,
                                         solver->_getTetSpecCount(oTet()->idx(), gidxion) + count);
            }
            solver->_setTetSpecCount(iTet()->idx(),
                                     gidxion,
                                     solver->_getTetSpecCount(iTet()->idx(), gidxion) - count);
        }
    }

    return current;
}

}  // namespace steps::tetode
