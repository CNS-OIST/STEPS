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
#include <iostream>
#include <fstream>
#include <sstream>
// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/mpi/tetopsplit/vdeptrans.hpp"
#include "steps/mpi/tetopsplit/tri.hpp"
#include "steps/mpi/tetopsplit/tet.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/mpi/tetopsplit/tetopsplit.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::VDepTrans::VDepTrans(ssolver::VDepTransdef * vdtdef, smtos::Tri * tri)
: KProc()
, pVDepTransdef(vdtdef)
, pTri(tri)
, localUpdVec()
, remoteUpdVec()
{
    assert (pVDepTransdef != 0);
    assert (pTri != 0);
    type = KP_VDEPTRANS;
}

////////////////////////////////////////////////////////////////////////////////

smtos::VDepTrans::~VDepTrans(void)
{
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::checkpoint(std::fstream & cp_file)
{
    cp_file.write((char*)&rExtent, sizeof(uint));
    cp_file.write((char*)&pFlags, sizeof(uint));

    cp_file.write((char*)&(crData.recorded), sizeof(bool));
    cp_file.write((char*)&(crData.pow), sizeof(int));
    cp_file.write((char*)&(crData.pos), sizeof(unsigned));
    cp_file.write((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::restore(std::fstream & cp_file)
{
    cp_file.read((char*)&rExtent, sizeof(uint));
    cp_file.read((char*)&pFlags, sizeof(uint));

    cp_file.read((char*)&(crData.recorded), sizeof(bool));
    cp_file.read((char*)&(crData.pow), sizeof(int));
    cp_file.read((char*)&(crData.pos), sizeof(unsigned));
    cp_file.read((char*)&(crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::reset(void)
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::setupDeps(void)
{
    assert(pTri->getInHost());
    std::set<smtos::KProc*> updset;

    uint nkprocs = pTri->countKProcs();
    uint startKProcIdx = pTri->getStartKProcIdx();
    
    // check if sk KProc in pTri depends on spec in pTri
    for (uint sk = 0; sk < nkprocs; sk++)
    {
        if (pTri->KProcDepSpecTri(sk, pTri, pVDepTransdef->srcchanstate()) == true)
            updset.insert(pTri->getKProc(sk));
        if (pTri->KProcDepSpecTri(sk, pTri, pVDepTransdef->dstchanstate()) == true)
            updset.insert(pTri->getKProc(sk));
    }

    localUpdVec.assign(updset.begin(), updset.end());

}

////////////////////////////////////////////////////////////////////////////////

bool smtos::VDepTrans::depSpecTet(uint gidx, smtos::WmVol * tet)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::VDepTrans::depSpecTri(uint gidx, smtos::Tri * triangle)
{
    if (triangle != pTri) return false;
    return (pVDepTransdef->dep(gidx) != ssolver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::VDepTrans::rate(steps::mpi::tetopsplit::TetOpSplitP * solver)
{
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint vdtlidx = pdef->vdeptransG2L(pVDepTransdef->gidx());
    // Fetch the local index of the srcchannel
    uint srclidx = pdef->vdeptrans_srcchanstate(vdtlidx);

    double n = static_cast<double>(pTri->pools()[srclidx]);
    double v = solver->getTriV(pTri->idx());
    double ra = pVDepTransdef->getVDepRate(v);

    return ra*n;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::apply(steps::rng::RNG * rng, double dt, double simtime, double period)
{
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint lidx = pdef->vdeptransG2L(pVDepTransdef->gidx());

    uint src = pdef->vdeptrans_srcchanstate(lidx);
    uint dst = pdef->vdeptrans_dstchanstate(lidx);


    uint nocs = pdef->countOhmicCurrs();
    for (uint oc = 0; oc < nocs; ++oc)
    {
        uint oc_cs = pdef->ohmiccurr_chanstate(oc);
        if (oc_cs == src)
        {
            if (pTri->clamped(src) == true) continue;
            pTri->setOCchange(oc, src, dt, simtime);
        }
        else if (oc_cs == dst)
        {
            if (pTri->clamped(dst) == true) continue;
            pTri->setOCchange(oc, dst, dt, simtime);
        }
    }

	if (pTri->clamped(src) == false)
	{
		uint nc = pTri->pools()[src];
		assert(nc >= 1);
		pTri->setCount(src,  (nc-1), period);
	}
	if (pTri->clamped(dst) == false)
	{
		uint nc = pTri->pools()[dst];
		assert(nc >= 0);
		pTri->setCount(dst,  (nc+1), period);
	}

    rExtent++;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<uint> const & smtos::VDepTrans::getRemoteUpdVec(int direction)
{
    return remoteUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::VDepTrans::getLocalUpdVec(int direction)
{
    return localUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::VDepTrans::resetOccupancies(void)
{
    
    pTri->resetPoolOccupancy();
}

////////////////////////////////////////////////////////////////////////////////

// END
