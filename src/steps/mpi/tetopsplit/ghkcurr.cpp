/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
#include <vector>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <limits>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/math/ghk.hpp"
#include "steps/mpi/tetopsplit/ghkcurr.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/wmvol.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;
namespace sm = steps::math;

////////////////////////////////////////////////////////////////////////////////

smtos::GHKcurr::GHKcurr(ssolver::GHKcurrdef * ghkdef, smtos::Tri * tri)
: KProc()
, pGHKcurrdef(ghkdef)
, pTri(tri)
, localUpdVec()
, remoteUpdVec()
, pEffFlux(true)
{
    assert (pGHKcurrdef != 0);
    assert (pTri != 0);
}

////////////////////////////////////////////////////////////////////////////////

smtos::GHKcurr::~GHKcurr(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));
    cp_file.write((char*)&pEffFlux, sizeof(bool));

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&rExtent, sizeof(uint));
    cp_file.read((char*)&pFlags, sizeof(uint));
    cp_file.read((char*)&pEffFlux, sizeof(bool));

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::reset(void)
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    setActive(true);
    pEffFlux = true;    //TODO: come back to this and check if rate needs to be recalculated here
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::setupDeps(void)
{
    assert(pTri->getInHost());
    std::set<smtos::KProc*> local;

    // The only concentration changes for a GHK current event are in the outer
    // and inner volume. The flux can involve movement of ion from either
    // compartment to the other- depnding on direction of flux

    WmVol * itet = pTri->iTet();
    WmVol * otet = pTri->oTet();
    assert (itet != 0);
    
    if (itet->getHost() != pTri->getHost()) {
        std::ostringstream os;
        os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
        throw steps::NotImplErr(os.str());
    }
    
    // The global species of the ion
    const uint gidxion = pGHKcurrdef->ion();

    // First check KProcs in the inner tetrahedron
    uint nkprocs = itet->countKProcs();
    uint startKProcIdx = itet->getStartKProcIdx();
    for (uint k = 0; k < nkprocs; k++)
    {
        if (itet->KProcDepSpecTet(k, itet, gidxion) == true) {
            local.insert(itet->getKProc(k));
        }

    }

    std::vector<smtos::Tri *>::const_iterator tri_end = itet->nexttriEnd();
    for (std::vector<smtos::Tri *>::const_iterator tri = itet->nexttriBegin();
             tri != tri_end; ++tri)
    {
        if ((*tri) == 0) continue;
        
        if (itet->getHost() != (*tri)->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << (*tri)->idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
            throw steps::NotImplErr(os.str());
        }

        nkprocs = (*tri)->countKProcs();
        startKProcIdx = (*tri)->getStartKProcIdx();;
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            if ((*tri)->KProcDepSpecTet(sk, itet, gidxion) == true)
                local.insert((*tri)->getKProc(sk));
        }
    }

    if (otet != 0)
    {
        if (otet->getHost() != pTri->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron " << otet->idx()  << " belong to different hosts.\n";
            throw steps::NotImplErr(os.str());
        }
        
        // Now check KProcs in the outer tetrahedron
        nkprocs = otet->countKProcs();
        startKProcIdx = otet->getStartKProcIdx();
        for (uint k = 0; k < nkprocs; k++)
        {
            if (otet->KProcDepSpecTet(k, otet, gidxion) == true) {
                local.insert(otet->getKProc(k));
            }
        }

        tri_end = otet->nexttriEnd();
        for (std::vector<smtos::Tri *>::const_iterator tri = otet->nexttriBegin();
            tri != tri_end; ++tri)
        {
            if ((*tri) == 0) continue;
            if (otet->getHost() != (*tri)->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << (*tri)->idx() << " and its compartment tetrahedron " << otet->idx()  << " belong to different hosts.\n";
                throw steps::NotImplErr(os.str());
            }
            nkprocs = (*tri)->countKProcs();
            startKProcIdx = (*tri)->getStartKProcIdx();;
            for (uint sk = 0; sk < nkprocs; sk++)
            {
                if ((*tri)->KProcDepSpecTet(sk, otet, gidxion) == true)
                    local.insert((*tri)->getKProc(sk));
            }
        }
    }

    localUpdVec.assign(local.begin(), local.end());
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::GHKcurr::depSpecTet(uint gidx, smtos::WmVol * tet)
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

bool smtos::GHKcurr::depSpecTri(uint gidx, smtos::Tri * triangle)
{
    if (triangle != pTri) return false;
    return (pGHKcurrdef->dep(gidx) != ssolver::DEP_NONE);

}

////////////////////////////////////////////////////////////////////////////////
/*
double smtos::GHKcurr::rate(double v, double T)
{
    // The rate comes from the GHK flux equation. The flux equation
    // returns a single-channel current which must be converted to a rate,
    // remembering that flux can be positive or negative (bi-directional)
    const uint gidxion = pGHKcurrdef->ion();

    double iconc = pTri->iTet()->conc(gidxion);
    double oconc = pTri->oTet()->conc(gidxion);
    double flux = sm::GHKcurrent(pGHKcurrdef->perm(), v, pGHKcurrdef->valence(),
                                 T, iconc, oconc);

    // Note: For a positive flux, this could be an efflux of +ve cations,
    // or an influx of -ve anions. Need to check the valence.

    // Get the rate of ion flux, remembering valence may by other than 1.
    double rt = flux/(sm::E_CHARGE * pGHKcurrdef->valence());
    // Now a positive rate is always an efflux and a negative rate is an influx
    // Set positive or negative flux flag TODO check the convention, could be reversed
    if (rt >= 0.0) pEffFlux = true;
    else pEffFlux = false;

    // Find the number of available channel states
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());
    // Fetch the local index of the channelstate
    uint cslidx = pdef->ghkcurr_chanstate(ghklidx);
    double n = static_cast<double>(pTri->pools()[cslidx]);

    return fabs(rt) * n;
}
*/

double smtos::GHKcurr::rate(steps::mpi::tetopsplit::TetOpSplitP * solver)
{
    // The rate comes from the GHK flux equation. The flux equation
    // returns a single-channel current which must be converted to a rate,
    // remembering that flux can be positive or negative (bi-directional)
    const uint gidxion = pGHKcurrdef->ion();
    double voconc = pGHKcurrdef->voconc();

    // Get concentrations in Molar units: convert to Mol/m^3
    double iconc = (pTri->iTet()->conc(gidxion))*1.0e3;
    double oconc = 0.0;

    if (voconc < 0.0)  oconc = (pTri->oTet()->conc(gidxion))*1.0e3;
    else  oconc = voconc*1.0e3;

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
    if (rt >= 0.0) setEffFlux(true);
    else setEffFlux(false);

    // Find the number of available channel states
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());
    // Fetch the local index of the channelstate
    uint cslidx = pdef->ghkcurr_chanstate(ghklidx);
    double n = static_cast<double>(pTri->pools()[cslidx]);

    return fabs(rt) * n;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::apply(steps::rng::RNG * rng, double dt, double simtime, double period) {
    smtos::WmVol * itet = pTri->iTet();
    smtos::WmVol * otet = pTri->oTet();

    ssolver::Compdef * innercdef = itet->compdef();
    ssolver::Compdef * outercdef = 0;
    if (otet) outercdef = otet->compdef();

    // Fetch the global index of the ion and the valence
    const uint gidxion = pGHKcurrdef->ion();
    const uint valence = pGHKcurrdef->valence();

    ssolver::Patchdef * pdef = pTri->patchdef();
    const uint ghklidx = pdef->ghkcurrG2L(pGHKcurrdef->gidx());

    assert (valence != 0);        // TODO: get rid of this when tested

    const bool realflux = pGHKcurrdef->realflux();
    double voconc = pGHKcurrdef->voconc();

    uint linneridx = innercdef->specG2L(gidxion);

    uint louteridx = std::numeric_limits<unsigned int>::max( );
    if (outercdef) louteridx = outercdef->specG2L(gidxion);

    if (efflux())
    {
        if (realflux)
        {
            if (itet->clamped(linneridx) == false) itet->incCount(linneridx,- 1, period);
            if (otet && voconc < 0.0)
            {
                if (otet->clamped(louteridx) == false) otet->incCount(louteridx, 1, period);
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, valence);
    }
    else
    {
        if (realflux)
        {
            if (itet->clamped(linneridx) == false) itet->incCount(linneridx, 1, period);
            if (otet && voconc < 0.0)
            {
                if (otet->clamped(louteridx) == false) otet->incCount(louteridx, -1, period);
            }
        }
        // A positive outward current is a positive current by convention
        pTri->incECharge(ghklidx, -valence );
    }

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::GHKcurr::getRemoteUpdVec(int direction)
{
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::GHKcurr::getLocalUpdVec(int direction)
{
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::GHKcurr::resetOccupancies(void)
{
    
    pTri->resetPoolOccupancy();
    
    // Update inner tet pools.
    smtos::WmVol * itet = pTri->iTet();
    if (itet != 0)
    {
        itet->resetPoolOccupancy();
    }
    
    smtos::WmVol * otet = pTri->oTet();
    if (otet != 0)
    {
        otet->resetPoolOccupancy();
    }
    
}

////////////////////////////////////////////////////////////////////////////////

// END
