/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/math/constants.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/vdeptrans.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::VDepTrans::VDepTrans(ssolver::VDepTransdef * vdtdef, stex::Tri * tri)
: 
 pVDepTransdef(vdtdef)
, pTri(tri)
, pUpdVec()
{
    AssertLog(pVDepTransdef != nullptr);
    AssertLog(pTri != nullptr);
}

////////////////////////////////////////////////////////////////////////////////

stex::VDepTrans::~VDepTrans()
= default;

////////////////////////////////////////////////////////////////////////////////

void stex::VDepTrans::checkpoint(std::fstream & cp_file)
{
    cp_file.write(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.write(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    cp_file.write(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.write(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.write(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.write(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::VDepTrans::restore(std::fstream & cp_file)
{
    cp_file.read(reinterpret_cast<char*>(&rExtent), sizeof(unsigned long long));
    cp_file.read(reinterpret_cast<char*>(&pFlags), sizeof(uint));

    cp_file.read(reinterpret_cast<char*>(&crData.recorded), sizeof(bool));
    cp_file.read(reinterpret_cast<char*>(&crData.pow), sizeof(int));
    cp_file.read(reinterpret_cast<char*>(&crData.pos), sizeof(unsigned));
    cp_file.read(reinterpret_cast<char*>(&crData.rate), sizeof(double));
}

////////////////////////////////////////////////////////////////////////////////

void stex::VDepTrans::reset()
{

    crData.recorded = false;
    crData.pow = 0;
    crData.pos = 0;
    crData.rate = 0.0;
    setActive(true);
}

////////////////////////////////////////////////////////////////////////////////

void stex::VDepTrans::setupDeps()
{
    std::set<stex::KProc*> updset;

    for (auto const& k : pTri->kprocs()) {
        if (k->depSpecTri(pVDepTransdef->srcchanstate(), pTri)) {
            updset.insert(k);
        } else if (k->depSpecTri(pVDepTransdef->dstchanstate(), pTri)) {
            updset.insert(k);
        }
    }

    pUpdVec.assign(updset.begin(), updset.end());

}

////////////////////////////////////////////////////////////////////////////////

bool stex::VDepTrans::depSpecTet(uint /*gidx*/, stex::WmVol * /*tet*/)
{
    return false;
}

////////////////////////////////////////////////////////////////////////////////

bool stex::VDepTrans::depSpecTri(uint gidx, stex::Tri * triangle)
{
    if (triangle != pTri) { return false;
}
    return (pVDepTransdef->dep(gidx) != ssolver::DEP_NONE);
}

////////////////////////////////////////////////////////////////////////////////

double stex::VDepTrans::rate(steps::tetexact::Tetexact * solver)
{
    ssolver::Patchdef * pdef = pTri->patchdef();
    uint vdtlidx = pdef->vdeptransG2L(pVDepTransdef->gidx());
    // Fetch the local index of the srcchannel
    uint srclidx = pdef->vdeptrans_srcchanstate(vdtlidx);

    auto n = static_cast<double>(pTri->pools()[srclidx]);
    double v = solver->getTriV(pTri->idx());
    double ra = pVDepTransdef->getVDepRate(v);

    return ra*n;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<stex::KProc*> const & stex::VDepTrans::apply(const rng::RNGptr &/*rng*/, double dt, double simtime)
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
            if (pTri->clamped(src)) continue;
            pTri->setOCchange(oc, src, dt, simtime);
        }
        else if (oc_cs == dst)
        {
            if (pTri->clamped(dst)) continue;
            pTri->setOCchange(oc, dst, dt, simtime);
        }
    }

    if (pTri->clamped(src) == false)
    {
        uint nc = pTri->pools()[src];
        AssertLog(nc >= 1);
        pTri->setCount(src,  (nc-1));
    }
    if (pTri->clamped(dst) == false)
    {
        uint nc = pTri->pools()[dst];
        pTri->setCount(dst,  (nc+1));
    }

    rExtent++;

    return pUpdVec;
}

////////////////////////////////////////////////////////////////////////////////

// END
