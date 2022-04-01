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
#include "tri.hpp"

#include <algorithm>
#include <cassert>
#include <cmath>
#include <fstream>
#include <functional>
#include <iostream>
#include <sstream>

#include <easylogging++.h>
#include <mpi.h>

#include "ghkcurr.hpp"
#include "sdiff.hpp"
#include "sreac.hpp"
#include "tet.hpp"
#include "tetopsplit.hpp"
#include "vdepsreac.hpp"
#include "vdeptrans.hpp"

#include "mpi/mpi_common.hpp"
#include "math/constants.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/sreacdef.hpp"
#include "util/error.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;
namespace ssolver = steps::solver;

////////////////////////////////////////////////////////////////////////////////

smtos::Tri::Tri(triangle_id_t idx, steps::solver::Patchdef *patchdef, double area,
                double l0, double l1, double l2, double d0, double d1, double d2,
                tetrahedron_id_t tetinner, tetrahedron_id_t tetouter,
                triangle_id_t tri0, triangle_id_t tri1, triangle_id_t tri2,
                int rank, int host_rank)
: pIdx(idx)
, pPatchdef(patchdef)
, pArea(area)
, pLengths()
, pDist()
, pTets()
, pNextTri()
, pECharge_last(nullptr)
, pECharge_accum(nullptr)
, pECharge_last_dt(0)
, pECharge_accum_dt(0)
, hostRank(host_rank)
, myRank(rank)
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

    pECharge_accum = new int[nghkcurrs];
    std::fill_n(pECharge_accum, nghkcurrs, 0);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    pOCchan_timeintg = new double[nohmcurrs];
    std::fill_n(pOCchan_timeintg, nohmcurrs, 0.0);
    pOCtime_upd = new double[nohmcurrs];
    std::fill_n(pOCtime_upd, nohmcurrs, 0.0);

    pPoolOccupancy = new double[nspecs];
    std::fill_n(pPoolOccupancy, nspecs, 0.0);
    pLastUpdate = new double[nspecs];
    std::fill_n(pLastUpdate, nspecs, 0.0);

    std::fill_n(pSDiffBndDirection, 3, false);
}

////////////////////////////////////////////////////////////////////////////////

smtos::Tri::~Tri()
{
    delete[] pPoolCount;
    delete[] pPoolFlags;
    delete[] pECharge;
    delete[] pECharge_last;
    delete[] pECharge_accum;
    delete[] pOCchan_timeintg;
    delete[] pOCtime_upd;
    delete[] pPoolOccupancy;
    delete[] pLastUpdate;

    KProcPVecCI e = pKProcs.end();
    for (std::vector<smtos::KProc *>::const_iterator i = pKProcs.begin();
         i != e; ++i) delete *i;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::checkpoint(std::fstream & cp_file)
{
    uint nspecs = patchdef()->countSpecs();
    cp_file.write(reinterpret_cast<char*>(pPoolCount), sizeof(uint) * nspecs);
    cp_file.write(reinterpret_cast<char*>(pPoolFlags), sizeof(uint) * nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    cp_file.write(reinterpret_cast<char*>(pECharge), sizeof(int) * nghkcurrs);
    cp_file.write(reinterpret_cast<char*>(pECharge_last), sizeof(int) * nghkcurrs);
    cp_file.write(reinterpret_cast<char*>(pECharge_accum), sizeof(int) * nghkcurrs);
    cp_file.write(reinterpret_cast<char*>(&pECharge_last_dt), sizeof(double));
    cp_file.write(reinterpret_cast<char*>(&pECharge_accum_dt), sizeof(double));

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    cp_file.write(reinterpret_cast<char*>(pOCchan_timeintg), sizeof(double) * nohmcurrs);
    cp_file.write(reinterpret_cast<char*>(pOCtime_upd), sizeof(double) * nohmcurrs);

    cp_file.write(reinterpret_cast<char*>(pSDiffBndDirection), sizeof(bool) * 3);

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::restore(std::fstream & cp_file)
{
    uint nspecs = patchdef()->countSpecs();
    cp_file.read(reinterpret_cast<char*>(pPoolCount), sizeof(uint) * nspecs);
    cp_file.read(reinterpret_cast<char*>(pPoolFlags), sizeof(uint) * nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    cp_file.read(reinterpret_cast<char*>(pECharge), sizeof(int) * nghkcurrs);
    cp_file.read(reinterpret_cast<char*>(pECharge_last), sizeof(int) * nghkcurrs);
    cp_file.read(reinterpret_cast<char*>(pECharge_accum), sizeof(int) * nghkcurrs);
    cp_file.read(reinterpret_cast<char*>(&pECharge_last_dt), sizeof(double));
    cp_file.read(reinterpret_cast<char*>(&pECharge_accum_dt), sizeof(double));

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    cp_file.read(reinterpret_cast<char*>(pOCchan_timeintg), sizeof(double) * nohmcurrs);
    cp_file.read(reinterpret_cast<char*>(pOCtime_upd), sizeof(double) * nohmcurrs);

    cp_file.read(reinterpret_cast<char*>(pSDiffBndDirection), sizeof(bool) * 3);

}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setInnerTet(smtos::WmVol * t)
{
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setOuterTet(smtos::WmVol * t)
{
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setSDiffBndDirection(uint i)
{
    AssertLog(i < 3);

    pSDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setNextTri(uint i, smtos::Tri * t)
{
    AssertLog(i <= 2);

    pNextTri[i]= t;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setupKProcs(smtos::TetOpSplitP * tex, bool efield)
{
    hasEfield = efield;
    startKProcIdx = tex->countKProcs();
    uint j = 0;

    nKProcs = pPatchdef->countSReacs()+pPatchdef->countSurfDiffs();
    if (hasEfield) {
        nKProcs += (pPatchdef->countVDepTrans() + pPatchdef->countVDepSReacs() + pPatchdef->countGHKcurrs());
}

    if (hostRank == myRank) {
        pKProcs.resize(nKProcs);


        // Create surface reaction kprocs
        uint nsreacs = patchdef()->countSReacs();
        for (uint i=0; i < nsreacs; ++i)
        {
            ssolver::SReacdef * srdef = patchdef()->sreacdef(i);
            auto * sr = new SReac(srdef, this);
            AssertLog(sr != nullptr);
            pKProcs[j++] = sr;
            uint idx = tex->addKProc(sr);
            sr->setSchedIDX(idx);
        }

        uint nsdiffs = patchdef()->countSurfDiffs();
        for (uint i=0; i < nsdiffs; ++i)
        {
            ssolver::Diffdef * sddef = patchdef()->surfdiffdef(i);
            auto * sd = new SDiff(sddef, this);
            AssertLog(sd != nullptr);
            pKProcs[j++] = sd;
            uint idx = tex->addKProc(sd);
            sd->setSchedIDX(idx);
            tex->addSDiff(sd);
        }

        // REMINDER: hasEfield not ready
        if (hasEfield)
        {
            uint nvdtrans = patchdef()->countVDepTrans();
            for (uint i=0; i < nvdtrans; ++i)
            {
                ssolver::VDepTransdef * vdtdef = patchdef()->vdeptransdef(i);
                auto * vdt = new VDepTrans(vdtdef, this);
                AssertLog(vdt != nullptr);
                pKProcs[j++] = vdt;
                uint idx = tex->addKProc(vdt);
                vdt->setSchedIDX(idx);
            }

            uint nvdsreacs = patchdef()->countVDepSReacs();
            for (uint i=0; i < nvdsreacs; ++i)
            {
                ssolver::VDepSReacdef * vdsrdef = patchdef()->vdepsreacdef(i);
                auto * vdsr = new VDepSReac(vdsrdef, this);
                AssertLog(vdsr != nullptr);
                pKProcs[j++] = vdsr;
                uint idx = tex->addKProc(vdsr);
                vdsr->setSchedIDX(idx);
            }

            uint nghkcurrs = patchdef()->countGHKcurrs();
            for (uint i=0; i < nghkcurrs; ++i)
            {
                ssolver::GHKcurrdef * ghkdef = patchdef()->ghkcurrdef(i);
                auto * ghk = new GHKcurr(ghkdef, this);
                AssertLog(ghk != nullptr);
                pKProcs[j++] = ghk;
                uint idx = tex->addKProc(ghk);
                ghk->setSchedIDX(idx);
            }
        }
    }
    else {
        pKProcs.resize(0);
        for (uint k = 0; k < nKProcs; k++) {
            tex->addKProc(nullptr);
        }
    }
}
////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setupDeps()
{
    if (myRank != hostRank) { return;
}
    for (auto& kp : pKProcs) {
        kp->setupDeps();
    }

    bool has_remote_neighbors = false;
    uint nspecs = patchdef()->countSpecs();
    for (uint i = 0; i < 3; ++i)
    {
        // Fetch next triangles, if it exists.
        smtos::Tri * next = nextTri(i);
        if (next == nullptr) { continue;
}
        if (next->getHost() != getHost()) {
            has_remote_neighbors = true;
            break;
        }
    }
    if (!has_remote_neighbors) {
        localSpecUpdKProcs.clear();
        return;
    }
    localSpecUpdKProcs.resize(nspecs);

    for (uint slidx = 0; slidx < nspecs; slidx++) {
        uint sgidx = patchdef()->specL2G(slidx);
        // check dependency of kprocs in the same tri
        uint nkprocs = countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++)
        {
            // Check locally.
            if (KProcDepSpecTri(sk, this, sgidx)) {
                localSpecUpdKProcs[slidx].push_back(getKProc(sk));
            }
        }

        // Check the neighbouring tetrahedrons.
        smtos::WmVol * itet = iTet();
        if (itet != nullptr)
        {
            if (getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << idx() << " and its compartment tetrahedron " << itet->idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = itet->countKProcs();
            for (uint k = 0; k < nkprocs; k++)
            {
                if (itet->KProcDepSpecTri(k, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(itet->getKProc(k));
                }
            }
        }

        smtos::WmVol * otet = oTet();
        if (otet != nullptr)
        {
            if (getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << idx() << " and its compartment tetrahedron " << otet->idx()  << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = otet->countKProcs();
            for (uint k = 0; k < nkprocs; k++)
            {
                if (otet->KProcDepSpecTri(k, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(otet->getKProc(k));
                }
            }
        }

    }
}

////////////////////////////////////////////////////////////////////////////////

bool smtos::Tri::KProcDepSpecTet(uint kp_lidx, smtos::WmVol* kp_container,  uint spec_gidx)
{
    // if kp is surf reaction
    uint remain = kp_lidx;
    if (remain < pPatchdef->countSReacs()) {
        ssolver::SReacdef * srdef = patchdef()->sreacdef(remain);
        if (kp_container == iTet())
        {
            return (srdef->dep_I(spec_gidx) != ssolver::DEP_NONE);
        }
        else if (kp_container == oTet())
        {
            return (srdef->dep_O(spec_gidx) != ssolver::DEP_NONE);
        }
        return false;
    }
    remain -= pPatchdef->countSReacs();
    // if kp is surface diff
    if (remain < pPatchdef->countSurfDiffs()) {
        return false;
    }
    remain -= pPatchdef->countSurfDiffs();

    // REMINDER: hasEfield not ready
    if (hasEfield)
    {
        // VDepTrans
        if (remain < pPatchdef->countVDepTrans()) {
            return false;
        }
        remain -= pPatchdef->countVDepTrans();

        // VDepSReac
        if (remain < pPatchdef->countVDepSReacs()) {
            ssolver::VDepSReacdef * vdsrdef = patchdef()->vdepsreacdef(remain);
            if (kp_container == iTet())
            {
                return (vdsrdef->dep_I(spec_gidx) != ssolver::DEP_NONE);
            }
            else if (kp_container == oTet())
            {
                return (vdsrdef->dep_O(spec_gidx) != ssolver::DEP_NONE);
            }
            return false;
        }
        remain -= pPatchdef->countVDepSReacs();

        // GHK
        if (remain < pPatchdef->countGHKcurrs()) {
            ssolver::GHKcurrdef * ghkdef = patchdef()->ghkcurrdef(remain);
            if (kp_container == iTet() )
            {
                return (ghkdef->dep_v(spec_gidx) != ssolver::DEP_NONE);
            }
            else if (kp_container == oTet())
            {
                if (ghkdef->voconc() < 0.0)
                {
                    return (ghkdef->dep_v(spec_gidx) != ssolver::DEP_NONE);
                }
                else
                {
                    return false;
                }
            }

            return false;
        }
    }

    AssertLog(false);
}


////////////////////////////////////////////////////////////////////////////////

bool smtos::Tri::KProcDepSpecTri(uint kp_lidx, smtos::Tri* kp_container, uint spec_gidx)
{
    // if kp is surf reaction
    uint remain = kp_lidx;
    if (remain < pPatchdef->countSReacs()) {
        if (kp_container != this) { return false;
}
        ssolver::SReacdef * srdef = patchdef()->sreacdef(remain);
        return (srdef->dep_S(spec_gidx) != ssolver::DEP_NONE);
    }
    remain -= pPatchdef->countSReacs();

    // if kp is surface diff
    if (remain < pPatchdef->countSurfDiffs()) {
        if (kp_container != this) { return false;
    }
        ssolver::Diffdef * sddef = patchdef()->surfdiffdef(remain);
        return spec_gidx == sddef->lig();
    }

    remain -= pPatchdef->countSurfDiffs();

    // REMINDER: hasEfield not ready
    if (hasEfield)
    {
        // VDepTrans
        if (remain < pPatchdef->countVDepTrans()) {
            ssolver::VDepTransdef * vdtdef = patchdef()->vdeptransdef(remain);
            if (kp_container != this) { return false;
}
            return (vdtdef->dep(spec_gidx) != ssolver::DEP_NONE);
        }
        remain -= pPatchdef->countVDepTrans();

        // VDepSReac
        if (remain < pPatchdef->countVDepSReacs()) {
            ssolver::VDepSReacdef * vdsrdef = patchdef()->vdepsreacdef(remain);
            if (kp_container != this) { return false;
}
            return (vdsrdef->dep_S(spec_gidx) != ssolver::DEP_NONE);
        }
        remain -= pPatchdef->countVDepSReacs();

        // GHK
        if (remain < pPatchdef->countGHKcurrs()) {
            ssolver::GHKcurrdef * ghkdef = patchdef()->ghkcurrdef(remain);
            if (kp_container != this) { return false;
}
            return (ghkdef->dep(spec_gidx) != ssolver::DEP_NONE);
        }
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::reset()
{
    uint nspecs = patchdef()->countSpecs();
    std::fill_n(pPoolCount, nspecs, 0);
    std::fill_n(pPoolFlags, nspecs, 0);

    for (auto kproc: pKProcs) {
        kproc->reset();
    }


    uint nghkcurrs = pPatchdef->countGHKcurrs();
    std::fill_n(pECharge, nghkcurrs, 0);
    std::fill_n(pECharge_last, nghkcurrs, 0);
    std::fill_n(pECharge_accum, nghkcurrs, 0);
    pECharge_last_dt = 0;
    pECharge_accum_dt = 0;

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    std::fill_n(pOCchan_timeintg, nohmcurrs, 0.0);
    std::fill_n(pOCtime_upd, nohmcurrs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::resetECharge(double dt, double efdt, double t)
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();

    for (uint i=0; i < nghkcurrs; ++i) {
        pECharge_accum[i] += pECharge[i];
    }
    pECharge_accum_dt += dt;

    if (pECharge_accum_dt >= efdt or
        (efdt - pECharge_accum_dt) <=
            std::numeric_limits<double>::epsilon() * t * 8) {
        // Swap arrays
        std::swap(pECharge_last, pECharge_accum);

        // reset accumulation array and dt values
        std::fill_n(pECharge_accum, nghkcurrs, 0);

        pECharge_last_dt = pECharge_accum_dt;
        pECharge_accum_dt = 0;
    }

    std::fill_n(pECharge, nghkcurrs, 0);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::resetOCintegrals()
{
    uint nocs = patchdef()->countOhmicCurrs();
    std::fill_n(pOCchan_timeintg, nocs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::incECharge(uint lidx, int charge)
{
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);
    pECharge[lidx]+=charge;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setCount(uint lidx, uint count, double period)
{
    // tri count doesn't care about global
    AssertLog(lidx < patchdef()->countSpecs());
	double oldcount = pPoolCount[lidx];
	pPoolCount[lidx] = count;
	if (period == 0.0) { return;
}
	// Count has changed,
	double lastupdate = pLastUpdate[lidx];
	AssertLog(period >= lastupdate);
	pPoolOccupancy[lidx] += oldcount*(period-lastupdate);

	pLastUpdate[lidx] = period;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::incCount(uint lidx, int inc, double period, bool local_change)
{
    // can change count locally or remotely, if change locally, check if the change require sync,
    // if change remotely, register the change in pSol, sync of remote change will be register when
    // the change is applied in remote host via this function
    AssertLog(lidx < patchdef()->countSpecs());

	// count changed by diffusion
    if (hostRank != myRank && !local_change)
    {
        if (inc <= 0) {
            std::ostringstream os;
            os << "Try to change molecule " << lidx << " by " << inc << "\n";
            os << "Fail because molecule change of receiving end should always be non-negative.\n";
            ProgErrLog(os.str());
        }
        bufferLocations[lidx] = pSol->registerRemoteMoleculeChange(hostRank, bufferLocations[lidx], SUB_TRI, pIdx.get(), lidx, inc);
    }
    // local change by reac or diff
    else {
        // need to check if the change costs negative molecule count
		double oldcount = pPoolCount[lidx];
		AssertLog(oldcount + inc >= 0.0);
		pPoolCount[lidx] += inc;

		if (period == 0.0 || local_change) { return;
}
		// Count has changed,
		double lastupdate = pLastUpdate[lidx];
		AssertLog(period >= lastupdate);
		pPoolOccupancy[lidx] += oldcount*(period-lastupdate);

		pLastUpdate[lidx] = period;

    }
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setOCchange(uint oclidx, uint slidx, double dt, double simtime)
{
    // NOTE: simtime is BEFORE the update has taken place

    AssertLog(oclidx < patchdef()->countOhmicCurrs());
    AssertLog(slidx < patchdef()->countSpecs());

    // A channel state relating to an ohmic current has changed it's
    // number.
    double integral = pPoolCount[slidx]*((simtime+dt) - pOCtime_upd[oclidx]);

    AssertLog(integral >= 0.0);

    pOCchan_timeintg[oclidx] += integral;
    pOCtime_upd[oclidx] = simtime+dt;
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setClamped(uint lidx, bool clamp)
{
    if (clamp) { pPoolFlags[lidx] |= CLAMPED;
    } else { pPoolFlags[lidx] &= ~CLAMPED;
}
}

////////////////////////////////////////////////////////////////////////////////

smtos::SReac * smtos::Tri::sreac(uint lidx) const
{
    AssertLog(lidx < patchdef()->countSReacs());
    return dynamic_cast<smtos::SReac*>(pKProcs[lidx]);
}

////////////////////////////////////////////////////////////////////////////////
int smtos::Tri::getTriDirection(triangle_id_t tidx)
{
    for (uint i = 0; i < 3; i++) {
        if (pTris[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

smtos::SDiff * smtos::Tri::sdiff(uint lidx) const
{
    AssertLog(lidx < patchdef()->countSurfDiffs());
    return dynamic_cast<smtos::SDiff*>(pKProcs[patchdef()->countSReacs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

smtos::VDepTrans * smtos::Tri::vdeptrans(uint lidx) const
{
    AssertLog(lidx < patchdef()->countVDepTrans());
    return dynamic_cast<smtos::VDepTrans*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + lidx]);
}

////////////////////////////////////////////////////////////////////////////////

smtos::VDepSReac * smtos::Tri::vdepsreac(uint lidx) const
{
    AssertLog(lidx < patchdef()->countVDepSReacs());
    return dynamic_cast<smtos::VDepSReac*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + patchdef()->countVDepTrans() + lidx]);

}

////////////////////////////////////////////////////////////////////////////////

smtos::GHKcurr * smtos::Tri::ghkcurr(uint lidx) const
{
    AssertLog(lidx < patchdef()->countGHKcurrs());
    return dynamic_cast<smtos::GHKcurr*>(pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + patchdef()->countVDepTrans() + patchdef()->countVDepSReacs() + lidx]);

}

////////////////////////////////////////////////////////////////////////////////

double smtos::Tri::getGHKI(uint lidx) const
{
    if (pECharge_last_dt == 0)
        return 0;

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);

    int efcharge = pECharge_last[lidx];
    auto efcharged = static_cast<double>(efcharge);

    return ((efcharged*steps::math::E_CHARGE)/pECharge_last_dt);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Tri::getGHKI() const
{
    if (pECharge_last_dt == 0)
        return 0;

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge=0;
    for (uint i =0; i < nghkcurrs; ++i)
    {
        efcharge += pECharge_last[i];
    }

    auto efcharged = static_cast<double>(efcharge);

    return ((efcharged*steps::math::E_CHARGE)/pECharge_last_dt);
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Tri::computeI(double v, double dt, double simtime, double efdt)
{
    /*
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (uint i = 0; i < nocs; ++i)
    {
        ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(i);
        // The next is ok because Patchdef returns local index
        uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(i)];
        current += (n*ocdef->getG())*(v-ocdef->getERev());
    }
    */

    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (uint i = 0; i < nocs; ++i)
    {
        ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(i);
        // First calculate the last little bit up to the simtime
        double integral = pPoolCount[patchdef()->ohmiccurr_chanstate(i)]*(simtime - pOCtime_upd[i]);
        AssertLog(integral >= 0.0);
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
    auto efcharged = static_cast<double>(efcharge);

    // Convert charge to coulombs and find mean current
    current += ((efcharged*steps::math::E_CHARGE)/dt);
    resetECharge(dt, efdt, simtime);
    resetOCintegrals();

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Tri::getOhmicI(double v,double /*dt*/) const
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

double smtos::Tri::getOhmicI(uint lidx, double v,double /*dt*/) const
{
    AssertLog(lidx < patchdef()->countOhmicCurrs());
    ssolver::OhmicCurrdef * ocdef = patchdef()->ohmiccurrdef(lidx);
    uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(lidx)];

    return (n*ocdef->getG())*(v-ocdef->getERev());
}

////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool smtos::Tri::getInHost()
{
    return hostRank == myRank;
}

////////////////////////////////////////////////////////////////////////////////
void smtos::Tri::setHost(int host, int rank)
{
    hostRank = host;
    myRank = rank;
}
////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setSolver(steps::mpi::tetopsplit::TetOpSplitP* sol)
{
    pSol = sol;
}

////////////////////////////////////////////////////////////////////////////////

smtos::TetOpSplitP* smtos::Tri::solver()
{
    return pSol;
}
////////////////////////////////////////////////////////////////////////////////
double smtos::Tri::getPoolOccupancy(uint lidx)
{
	AssertLog(lidx < patchdef()->countSpecs());

	return pPoolOccupancy[lidx];
}

////////////////////////////////////////////////////////////////////////////////

double smtos::Tri::getLastUpdate(uint lidx)
{
	AssertLog(lidx < patchdef()->countSpecs());

	return pLastUpdate[lidx];
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::resetPoolOccupancy()
{
	uint nspecs = pPatchdef->countSpecs();
    std::fill_n(pPoolOccupancy, nspecs, 0.0);
    std::fill_n(pLastUpdate, nspecs, 0.0);
}

////////////////////////////////////////////////////////////////////////////////

std::vector<smtos::KProc*> const & smtos::Tri::getSpecUpdKProcs(uint slidx)
{
    return localSpecUpdKProcs[slidx];
}


////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::repartition(smtos::TetOpSplitP * tex, int rank, int host_rank)
{
    myRank = rank;
    hostRank = host_rank;

    // Delete reaction rules.
    KProcPVecCI e = pKProcs.end();
    for (KProcPVecCI i = pKProcs.begin(); i != e; ++i) delete *i;

    setupKProcs(tex);

    localSpecUpdKProcs.clear();
    bufferLocations.clear();
}

////////////////////////////////////////////////////////////////////////////////

void smtos::Tri::setupBufferLocations()
{
    uint nspecs = pPatchdef->countSpecs();
    bufferLocations.assign(nspecs, std::numeric_limits<uint>::max());
}


////////////////////////////////////////////////////////////////////////////////
//END
