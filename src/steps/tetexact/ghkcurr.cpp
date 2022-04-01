/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


// Standard library & STL headers.
#include <cmath>
#include <iostream>
#include <limits>
#include <vector>

// STEPS headers.
#include "ghkcurr.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "tri.hpp"
#include "wmvol.hpp"
#include "math/ghk.hpp"

#include <easylogging++.h>
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;
namespace sm = steps::math;

////////////////////////////////////////////////////////////////////////////////

stex::GHKcurr::GHKcurr(ssolver::GHKcurrdef * ghkdef, stex::Tri * tri)
:
 pGHKcurrdef(ghkdef)
, pTri(tri)
, pUpdVec()
, pEffFlux(true)
{
    AssertLog(pGHKcurrdef != nullptr);
    AssertLog(pTri != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

stex::GHKcurr::~GHKcurr() = default;

////////////////////////////////////////////////////////////////////////////////

void stex::GHKcurr::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));
    cp_file.write(reinterpret_cast<char*>(&pEffFlux), sizeof(bool));

    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::GHKcurr::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.read(reinterpret_cast<char*>(&pFlags), sizeof(uint));
    cp_file.read(reinterpret_cast<char*>(&pEffFlux), sizeof(bool));

    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::GHKcurr::reset()
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    setActive(true);
    pEffFlux = true;    //TODO: come back to this and check if rate needs to be recalculated here
}

////////////////////////////////////////////////////////////////////////////////

void stex::GHKcurr::setupDeps()
{
    std::set<stex::KProc*> updset;

    // The only concentration changes for a GHK current event are in the outer
    // and inner volume. The flux can involve movement of ion from either
    // compartment to the other- depnding on direction of flux

    WmVol * itet = pTri->iTet();
    WmVol * otet = pTri->oTet();
    AssertLog(itet != nullptr);
    // The global species of the ion
    const uint gidxion = pGHKcurrdef->ion();

    // First check KProcs in the inner tetrahedron
    for (auto const& k: itet->kprocs()) {
        if (k->depSpecTet(gidxion, itet))
            updset.insert(k);
    }

    for (auto const& tri : itet->nexttris()) {
        if (tri == nullptr) continue;

        for (auto const& k: tri->kprocs()) {
            if (k->depSpecTet(gidxion, itet))
                updset.insert(k);
        }
    }

    if (otet != nullptr)
    {
        // Now check KProcs in the outer tetrahedron
        for (auto const& k: otet->kprocs()) {
            if (k->depSpecTet(gidxion, otet))
                updset.insert(k);
        }

        for (auto const& tri : otet->nexttris()) {
            if (tri == nullptr) continue;

            for (auto const& k : tri->kprocs()) {
                if (k->depSpecTet(gidxion, otet))
                    updset.insert(k);
            }
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());
}

////////////////////////////////////////////////////////////////////////////////

bool stex::GHKcurr::depSpecTet(uint gidx, stex::WmVol * tet)
{
    // Two things here. Any changes of the count of the ion in the inner or
    // outer tetrahedron affect the single-channel rate. This will be dealt with
    // together for now in rate()
    // NOTE: DEP_NONE is defined in steps/sim/shared/types.hpp
    if (tet == pTri->iTet() )
    {
        return (pGHKcurrdef->dep_v(gidx) != ssolver::DEP_NONE);
    }
    else if (tet == pTri->oTet())
    {
        if (pGHKcurrdef->voconc() < 0.0)
        {
            return (pGHKcurrdef->dep_v(gidx) != ssolver::DEP_NONE);
        }
        else
        {
            return false;
        }
    }

    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::GHKcurr::depSpecTri(uint gidx, stex::Tri * triangle)
{
    if (triangle != pTri) { return false;
}
    return (pGHKcurrdef->dep(gidx) != ssolver::DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////

double stex::GHKcurr::rate(steps::tetexact::Tetexact * solver)
{
    // The rate comes from the GHK flux equation. The flux equation
    // returns a single-channel current which must be converted to a rate,
    // remembering that flux can be positive or negative (bi-directional)
    const uint gidxion = pGHKcurrdef->ion();
    double voconc = pGHKcurrdef->voconc();
    double oconc = 0.0;

    // Get concentrations in Molar units: convert to Mol/m^3
    double iconc = (pTri->iTet()->conc(gidxion))*1.0e3;

    if (voconc < 0.0) {  oconc = (pTri->oTet()->conc(gidxion))*1.0e3;
    } else {  oconc = voconc*1.0e3;
}

    double v = solver->getTriV(pTri->idx());
    double T = solver->getTemp();

    double flux = sm::GHKcurrent(pGHKcurrdef->perm(), v+pGHKcurrdef->vshift(), pGHKcurrdef->valence(),
                                 T, iconc, oconc);

    // Note: For a positive flux, this could be an efflux of +ve cations,
    // or an influx of -ve anions. Need to check the valence.

    // Get the rate of ion flux, remembering valence may by other than 1.
    double rt = flux/(sm::E_CHARGE * static_cast<double>(pGHKcurrdef->valence()));
    // Now a positive rate is always an efflux and a negative rate is an influx
    // Set positive or negative flux flag
    setEffFlux(rt >= 0.0);

    // Find the number of available channel states
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());
    // Fetch the local index of the channelstate
    uint cslidx = pdef->ghkcurr_chanstate(ghklidx);
    auto n = static_cast<double>(pTri->pools()[cslidx]);

    return fabs(rt) * n;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::GHKcurr::apply(const rng::RNGptr &/*rng*/, double /*dt*/, double /*simtime*/)
{
    stex::WmVol * itet = pTri->iTet();
    stex::WmVol * otet = pTri->oTet();

    ssolver::Compdef * innercdef = itet->compdef();
    ssolver::Compdef * outercdef = nullptr;
    if (otet) outercdef = otet->compdef();

    // Fetch the global index of the ion and the valence
    const auto gidxion = pGHKcurrdef->ion();
    const auto valence = pGHKcurrdef->valence();

    ssolver::Patchdef * pdef = pTri->patchdef();
    const uint ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());

    AssertLog(valence != 0);        // TODO: get rid of this when tested

    const bool realflux = pGHKcurrdef->realflux();
    double voconc = pGHKcurrdef->voconc();

    uint linneridx = innercdef->specG2L(gidxion);

    uint louteridx = std::numeric_limits<unsigned int>::max( );
    if (outercdef) louteridx = outercdef->specG2L(gidxion);

    if (efflux())
    {
        if (realflux)
        {
            if (!itet->clamped(linneridx)) itet->incCount(linneridx, - 1);
            if (otet && voconc < 0.0)
            {
                if (!otet->clamped(louteridx)) otet->incCount(louteridx, 1);
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, valence);
    }
    else
    {
        if (realflux)
        {
            if (!itet->clamped(linneridx)) itet->incCount(linneridx, 1);
            if (otet && voconc < 0.0)
            {
                if (!otet->clamped(louteridx)) otet->incCount(louteridx, -1);
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, -valence );
    }

    rExtent++;

    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
