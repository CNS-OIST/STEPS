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
#include <cassert>
#include <cmath>
#include <algorithm>
#include <functional>
#include <iostream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/sreacdef.hpp"
#include "steps/solver/ohmiccurrdef.hpp"
#include "steps/tetexact/vdeptrans.hpp"
#include "steps/tetexact/vdepsreac.hpp"
#include "steps/tetexact/ghkcurr.hpp"
#include "steps/tetexact/sreac.hpp"
#include "steps/tetexact/sdiff.hpp"
#include "steps/tetexact/tet.hpp"
#include "steps/tetexact/tri.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/tetexact/tetexact.hpp"
#include "steps/math/constants.hpp"

#include "third_party/easylogging++.h"

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

stex::Tri::Tri(uint idx, steps::solver::Patchdef * patchdef, double area,
                double l0, double l1, double l2, double d0, double d1, double d2,
                 int tetinner, int tetouter, int tri0, int tri1, int tri2)
: pIdx(idx)
, pPatchdef(patchdef)
, pArea(area)
, pLengths()
, pDist()
, pInnerTet(0)
, pOuterTet(0)
, pTets()
, pNextTri()
, pPoolCount(0)
, pPoolFlags(0)
, pKProcs()
, pECharge(0)
, pECharge_last(0)
, pOCchan_timeintg(0)
, pOCtime_upd(0)
{
    assert(pPatchdef != 0);
    assert (pArea > 0.0);

    assert (l0 > 0.0 && l1 > 0.0 && l2 > 0.0);
    assert (d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0);

    pTets[0] = tetinner;
    pTets[1] = tetouter;

    pTris[0] = tri0;
    pTris[1] = tri1;
    pTris[2] = tri2;

    pNextTri[0] = 0;
    pNextTri[1] = 0;
    pNextTri[2] = 0;

    pLengths[0] = l0;
    pLengths[1] = l1;
    pLengths[2] = l2;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;

    uint nspecs = pPatchdef->countSpecs();
    pPoolCount = new uint[nspecs];
    pPoolFlags = new uint[nspecs];
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    pECharge = new int[nghkcurrs];
    std::fill_n(pECharge, nghkcurrs, 0);

    pECharge_last = new int[nghkcurrs];
    std::fill_n(pECharge_last, nghkcurrs, 0);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    pOCchan_timeintg = new double[nohmcurrs];
    std::fill_n(pOCchan_timeintg, nohmcurrs, 0.0);
    pOCtime_upd = new double[nohmcurrs];
    std::fill_n(pOCtime_upd, nohmcurrs, 0.0);

    std::fill_n(pDiffBndDirection, 3, false);
}

////////////////////////////////////////////////////////////////////////////////

stex::Tri::~Tri(void)
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
    delete[] pECharge;
    delete[] pOCchan_timeintg;
    delete[] pOCtime_upd;

    KProcPVecCI e = pKProcs.end();
    for (std::vector<stex::KProc *>::const_iterator i = pKProcs.begin();
         i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::checkpoint(std::fstream & cp_file)
{
    uint nspecs = patchdef()->countSpecs();
    cp_file.write((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.write((char*)pPoolFlags, sizeof(uint) * nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    cp_file.write((char*)pECharge, sizeof(int) * nghkcurrs);
    cp_file.write((char*)pECharge_last, sizeof(int) * nghkcurrs);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    cp_file.write((char*)pOCchan_timeintg, sizeof(double) * nohmcurrs);
    cp_file.write((char*)pOCtime_upd, sizeof(double) * nohmcurrs);

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::restore(std::fstream & cp_file)
{
    uint nspecs = patchdef()->countSpecs();
    cp_file.read((char*)pPoolCount, sizeof(uint) * nspecs);
    cp_file.read((char*)pPoolFlags, sizeof(uint) * nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    cp_file.read((char*)pECharge, sizeof(int) * nghkcurrs);
    cp_file.read((char*)pECharge_last, sizeof(int) * nghkcurrs);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    cp_file.read((char*)pOCchan_timeintg, sizeof(double) * nohmcurrs);
    cp_file.read((char*)pOCtime_upd, sizeof(double) * nohmcurrs);

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setInnerTet(stex::WmVol * t)
{
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setOuterTet(stex::WmVol * t)
{
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setDiffBndDirection(uint i)
{
    assert(i < 3);

    pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setNextTri(uint i, stex::Tri * t)
{
    assert (i <= 2);

    pNextTri[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setupKProcs(stex::Tetexact * tex, bool efield)
{
    uint kprocvecsize = pPatchdef->countSReacs()+pPatchdef->countSurfDiffs();
    if (efield == true)
        kprocvecsize += (pPatchdef->countVDepTrans() + pPatchdef->countVDepSReacs() + pPatchdef->countGHKcurrs());
    pKProcs.resize(kprocvecsize);

    uint j = 0;
    // Create surface reaction kprocs
    uint nsreacs = patchdef()->countSReacs();
    for (uint i=0; i < nsreacs; ++i)
    {
        ssolver::SReacdef * srdef = patchdef()->sreacdef(i);
        stex::SReac * sr = new SReac(srdef, this);
        assert(sr != 0);
        pKProcs[j++] = sr;
        tex->addKProc(sr);
    }

    uint nsdiffs = patchdef()->countSurfDiffs();
    for (uint i=0; i < nsdiffs; ++i)
    {
        ssolver::Diffdef * sddef = patchdef()->surfdiffdef(i);
        stex::SDiff * sd = new SDiff(sddef, this);
        assert(sd != 0);
        pKProcs[j++] = sd;
        tex->addKProc(sd);
    }


    if (efield == true)
    {
        uint nvdtrans = patchdef()->countVDepTrans();
        for (uint i=0; i < nvdtrans; ++i)
        {
            ssolver::VDepTransdef * vdtdef = patchdef()->vdeptransdef(i);
            stex::VDepTrans * vdt = new VDepTrans(vdtdef, this);
            assert(vdt != 0);
            pKProcs[j++] = vdt;
            tex->addKProc(vdt);
        }

        uint nvdsreacs = patchdef()->countVDepSReacs();
        for (uint i=0; i < nvdsreacs; ++i)
        {
            ssolver::VDepSReacdef * vdsrdef = patchdef()->vdepsreacdef(i);
            stex::VDepSReac * vdsr = new VDepSReac(vdsrdef, this);
            assert(vdsr != 0);
            pKProcs[j++] = vdsr;
            tex->addKProc(vdsr);
        }

        uint nghkcurrs = patchdef()->countGHKcurrs();
        for (uint i=0; i < nghkcurrs; ++i)
        {
            ssolver::GHKcurrdef * ghkdef = patchdef()->ghkcurrdef(i);
            stex::GHKcurr * ghk = new GHKcurr(ghkdef, this);
            assert(ghk != 0);
            pKProcs[j++] = ghk;
            tex->addKProc(ghk);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::reset(void)
{
    uint nspecs = patchdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);

    std::for_each(pKProcs.begin(), pKProcs.end(),
        std::mem_fun(&stex::KProc::reset));


    uint nghkcurrs = pPatchdef->countGHKcurrs();
    std::fill_n(pECharge, nghkcurrs, 0);
    std::fill_n(pECharge_last, nghkcurrs, 0);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    std::fill_n(pOCchan_timeintg, nohmcurrs, 0.0);
    std::fill_n(pOCtime_upd, nohmcurrs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::resetECharge(void)
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();

    for (uint i=0; i < nghkcurrs; ++i)
    {
        pECharge_last[i] = pECharge[i];
    }
    std::fill_n(pECharge, nghkcurrs, 0);

}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::resetOCintegrals(void)
{
    uint nocs = patchdef()->countOhmicCurrs();
    std::fill_n(pOCchan_timeintg, nocs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::incECharge(uint lidx, int charge)
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    assert(lidx < nghkcurrs);
    pECharge[lidx]+=charge;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setCount(uint lidx, uint count)
{
    assert (lidx < patchdef()->countSpecs());
    double oldcount = pPoolCount[lidx];
    double c = static_cast<double>(count);
    pPoolCount[lidx] = c;

    /* 16/01/10 IH: Counts no longer stored in patch object.
    // Now update the count in this tri's patch
    double diff = c - oldcount;
    double newcount = (patchdef()->pools()[lidx]) + diff;
    // Patchdef method will do the checking on the double argument (should be positive!)
    patchdef()->setCount(lidx, newcount);
    */
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::incCount(uint lidx, int inc)
{
    assert (lidx < patchdef()->countSpecs());
    pPoolCount[lidx] += inc;
    assert(pPoolCount[lidx] >= 0);
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setOCchange(uint oclidx, uint slidx, double dt, double simtime)
{
    // NOTE: simtime is BEFORE the update has taken place

    assert(oclidx < patchdef()->countOhmicCurrs());
    assert(slidx < patchdef()->countSpecs());

    // A channel state relating to an ohmic current has changed it's
    // number.
    double integral = pPoolCount[slidx]*((simtime+dt) - pOCtime_upd[oclidx]);
    assert(integral >= 0.0);

    pOCchan_timeintg[oclidx] += integral;
    pOCtime_upd[oclidx] = simtime+dt;
}

////////////////////////////////////////////////////////////////////////////////

void stex::Tri::setClamped(uint lidx, bool clamp)
{
    if (clamp == true) pPoolFlags[lidx] |= CLAMPED;
    else pPoolFlags[lidx] &= ~CLAMPED;
}

////////////////////////////////////////////////////////////////////////////////

stex::SReac * stex::Tri::sreac(uint lidx) const
{
    assert(lidx < patchdef()->countSReacs());
    return dynamic_cast<stex::SReac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////
int stex::Tri::getTriDirection(uint tidx)
{
    for (uint i = 0; i < 3; i++) {
        if (pTris[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

stex::SDiff * stex::Tri::sdiff(uint lidx) const
{
    assert(lidx < patchdef()->countSurfDiffs());
    return dynamic_cast<stex::SDiff*>(pKProcs[patchdef()->countSReacs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

stex::VDepTrans * stex::Tri::vdeptrans(uint lidx) const
{
    assert(lidx < patchdef()->countVDepTrans());
    return dynamic_cast<stex::VDepTrans*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

stex::VDepSReac * stex::Tri::vdepsreac(uint lidx) const
{
    assert(lidx < patchdef()->countVDepSReacs());
    return dynamic_cast<stex::VDepSReac*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + patchdef()->countVDepTrans() + lidx]);

}

////////////////////////////////////////////////////////////////////////////////

stex::GHKcurr * stex::Tri::ghkcurr(uint lidx) const
{
    assert(lidx < patchdef()->countGHKcurrs());
    return dynamic_cast<stex::GHKcurr*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + patchdef()->countVDepTrans() + patchdef()->countVDepSReacs() + lidx]);

}

////////////////////////////////////////////////////////////////////////////////

double stex::Tri::getGHKI(uint lidx, double dt) const
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    assert(lidx < nghkcurrs);

    int efcharge = pECharge_last[lidx];
    double efcharged = static_cast<double>(efcharge);

    return ((efcharged*steps::math::E_CHARGE)/dt);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tri::getGHKI(double dt) const
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge=0;
    for (uint i =0; i < nghkcurrs; ++i)
    {
        efcharge += pECharge_last[i];
    }

    double efcharged = static_cast<double>(efcharge);

    return ((efcharged*steps::math::E_CHARGE)/dt);
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tri::computeI(double v, double dt, double simtime)
{
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (uint i = 0; i < nocs; ++i)
    {
        ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(i);
        // First calculate the last little bit up to the simtime
        double integral = pPoolCount[patchdef()->ohmiccurr_chanstate(i)]*(simtime - pOCtime_upd[i]);
        assert(integral >= 0.0);
        pOCchan_timeintg[i] += integral;
        pOCtime_upd[i] = simtime;

        // Find the mean number of channels open over the dt
        double n = pOCchan_timeintg[i]/dt;
        current += (n*ocdef->getG())*(v-ocdef->getERev());
    }
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge=0;
    for (uint i =0; i < nghkcurrs; ++i)
    {
        efcharge += pECharge[i];
    }
    
    // The contribution from GHK charge movement.
    double efcharged = static_cast<double>(efcharge);

    // Convert charge to coulombs and find mean current
    current += ((efcharged*steps::math::E_CHARGE)/dt);
    resetECharge();
    resetOCintegrals();

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tri::getOhmicI(double v,double dt) const
{
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (uint i = 0; i < nocs; ++i)
    {
        ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(i);
        // The next is ok because Patchdef returns local index
        uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(i)];
        current += (n*ocdef->getG())*(v-ocdef->getERev());
    }

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double stex::Tri::getOhmicI(uint lidx, double v,double dt) const
{
    assert(lidx < patchdef()->countOhmicCurrs());
    ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(lidx);
    uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(lidx)];

    return (n*ocdef->getG())*(v-ocdef->getERev());
}

////////////////////////////////////////////////////////////////////////////////

//END
