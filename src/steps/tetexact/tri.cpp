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

#include "tri.hpp"

// STEPS headers.
#include "ghkcurr.hpp"
#include "math/constants.hpp"
#include "sdiff.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/sreacdef.hpp"
#include "sreac.hpp"
#include "tet.hpp"
#include "tetexact.hpp"
#include "vdepsreac.hpp"

// logging
#include "util/error.hpp"
#include <easylogging++.h>

#include "util/checkpointing.hpp"

namespace steps::tetexact {

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
    , pInnerTet(nullptr)
    , pOuterTet(nullptr)
    , pArea(area)
    , pLengths()
    , pDist()
    , pECharge_last_dt(0)
    , pECharge_accum_dt(0) {
    AssertLog(pPatchdef != nullptr);
    AssertLog(pArea > 0.0);

    AssertLog(l0 > 0.0 && l1 > 0.0 && l2 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0);

    pTets[0] = tetinner;
    pTets[1] = tetouter;

    pTris[0] = tri0;
    pTris[1] = tri1;
    pTris[2] = tri2;

    std::fill(pNextTri.begin(), pNextTri.end(), nullptr);

    pLengths[0] = l0;
    pLengths[1] = l1;
    pLengths[2] = l2;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;

    uint nspecs = pPatchdef->countSpecs();
    pPoolCount.container().resize(nspecs);
    pPoolFlags.container().resize(nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    pECharge.container().resize(nghkcurrs);
    pECharge_last.container().resize(nghkcurrs);
    pECharge_accum.container().resize(nghkcurrs);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    pOCchan_timeintg.container().resize(nohmcurrs);
    pOCtime_upd.container().resize(nohmcurrs);

    std::fill(pSDiffBndDirection.begin(), pSDiffBndDirection.end(), false);
}

////////////////////////////////////////////////////////////////////////////////

Tri::~Tri() {
    auto e = pKProcs.end();
    for (auto i = pKProcs.begin(); i != e; ++i) {
        delete *i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tri::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pECharge_accum);
    util::checkpoint(cp_file, pECharge_accum_dt);
    util::checkpoint(cp_file, pECharge_last);
    util::checkpoint(cp_file, pECharge_last_dt);
    util::checkpoint(cp_file, pOCtime_upd);
    util::checkpoint(cp_file, pERev);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolFlags);
    util::restore(cp_file, pPoolCount);
    util::restore(cp_file, pECharge_accum);
    util::restore(cp_file, pECharge_accum_dt);
    util::restore(cp_file, pECharge_last);
    util::restore(cp_file, pECharge_last_dt);
    util::restore(cp_file, pOCtime_upd);
    util::restore(cp_file, pERev);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setInnerTet(WmVol* t) {
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setOuterTet(WmVol* t) {
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setSDiffBndDirection(uint i) {
    AssertLog(i < 3);

    pSDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setNextTri(uint i, Tri* t) {
    pNextTri.at(i) = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupKProcs(Tetexact* tex, bool efield) {
    uint kprocvecsize = pPatchdef->countSReacs() + pPatchdef->countSurfDiffs();
    if (efield) {
        kprocvecsize += (pPatchdef->countVDepSReacs() + pPatchdef->countGHKcurrs());
    }
    pKProcs.resize(kprocvecsize);

    uint j = 0;
    // Create surface reaction kprocs
    uint nsreacs = patchdef()->countSReacs();
    for (auto i: solver::sreac_local_id::range(nsreacs)) {
        solver::SReacdef* srdef = patchdef()->sreacdef(i);
        auto* sr = new SReac(srdef, this);
        AssertLog(sr != nullptr);
        pKProcs[j++] = sr;
        tex->addKProc(sr);
    }

    uint nsdiffs = patchdef()->countSurfDiffs();
    for (auto i: solver::surfdiff_local_id::range(nsdiffs)) {
        solver::SurfDiffdef* sddef = patchdef()->surfdiffdef(i);
        auto* sd = new SDiff(sddef, this);
        AssertLog(sd != nullptr);
        pKProcs[j++] = sd;
        tex->addKProc(sd);
    }

    if (efield) {
        uint nvdsreacs = patchdef()->countVDepSReacs();
        for (auto i: solver::vdepsreac_local_id::range(nvdsreacs)) {
            solver::VDepSReacdef* vdsrdef = patchdef()->vdepsreacdef(i);
            auto* vdsr = new VDepSReac(vdsrdef, this);
            AssertLog(vdsr != nullptr);
            pKProcs[j++] = vdsr;
            tex->addKProc(vdsr, true);
        }

        uint nghkcurrs = patchdef()->countGHKcurrs();
        for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
            solver::GHKcurrdef* ghkdef = patchdef()->ghkcurrdef(i);
            auto* ghk = new GHKcurr(ghkdef, this);
            AssertLog(ghk != nullptr);
            pKProcs[j++] = ghk;
            tex->addKProc(ghk, true);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tri::reset() {
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);

    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }

    std::fill(pECharge.begin(), pECharge.end(), 0);
    std::fill(pECharge_last.begin(), pECharge_last.end(), 0);
    std::fill(pECharge_accum.begin(), pECharge_accum.end(), 0);
    pECharge_last_dt = 0;
    pECharge_accum_dt = 0;

    std::fill(pOCchan_timeintg.begin(), pOCchan_timeintg.end(), 0.0);
    std::fill(pOCtime_upd.begin(), pOCtime_upd.end(), 0.0);

    pERev.clear();
}

////////////////////////////////////////////////////////////////////////////////

void Tri::resetECharge(double dt, double efdt, double t) {
    uint nghkcurrs = pPatchdef->countGHKcurrs();

    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        pECharge_accum[i] += pECharge[i];
    }
    pECharge_accum_dt += dt;

    if (pECharge_accum_dt >= efdt or
        (efdt - pECharge_accum_dt) <= std::numeric_limits<double>::epsilon() * t * 8) {
        // Swap arrays
        std::swap(pECharge_last.container(), pECharge_accum.container());

        // reset accumulation array and dt values
        std::fill(pECharge_accum.begin(), pECharge_accum.end(), 0);

        pECharge_last_dt = pECharge_accum_dt;
        pECharge_accum_dt = 0;
    }

    std::fill(pECharge.begin(), pECharge.end(), 0);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::resetOCintegrals() {
    std::fill(pOCchan_timeintg.begin(), pOCchan_timeintg.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::incECharge(solver::ghkcurr_local_id lidx, int charge) {
    pECharge.at(lidx) += charge;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setCount(solver::spec_local_id lidx, uint count) {
    pPoolCount.at(lidx) = count;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::incCount(solver::spec_local_id lidx, int inc) {
    pPoolCount.at(lidx) += inc;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setOCchange(solver::ohmiccurr_local_id oclidx,
                      solver::spec_local_id slidx,
                      double dt,
                      double simtime) {
    // NOTE: simtime is BEFORE the update has taken place

    AssertLog(oclidx < patchdef()->countOhmicCurrs());
    AssertLog(slidx < patchdef()->countSpecs());

    // A channel state relating to an ohmic current has changed it's
    // number.
    double integral = pPoolCount[slidx] * ((simtime + dt) - pOCtime_upd[oclidx]);
    AssertLog(integral >= 0.0);

    pOCchan_timeintg[oclidx] += integral;
    pOCtime_upd[oclidx] = simtime + dt;
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

void Tri::setClamped(solver::spec_local_id lidx, bool clamp) {
    if (clamp) {
        pPoolFlags[lidx] |= CLAMPED;
    } else {
        pPoolFlags[lidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

SReac& Tri::sreac(solver::sreac_local_id lidx) const {
    AssertLog(lidx < patchdef()->countSReacs());
    return *dynamic_cast<SReac*>(pKProcs[lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////
int Tri::getTriDirection(triangle_global_id tidx) const noexcept {
    for (uint i = 0; i < 3; i++) {
        if (pTris[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

SDiff& Tri::sdiff(solver::surfdiff_local_id lidx) const {
    AssertLog(lidx < patchdef()->countSurfDiffs());
    return *dynamic_cast<SDiff*>(pKProcs[patchdef()->countSReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac& Tri::vdepsreac(solver::vdepsreac_local_id lidx) const {
    AssertLog(lidx < patchdef()->countVDepSReacs());
    return *dynamic_cast<VDepSReac*>(
        pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr& Tri::ghkcurr(solver::ghkcurr_local_id lidx) const {
    AssertLog(lidx < patchdef()->countGHKcurrs());
    return *dynamic_cast<GHKcurr*>(
        pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() +
                patchdef()->countVDepSReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getGHKI(solver::ghkcurr_local_id lidx) const {
    if (pECharge_last_dt == 0) {
        return 0;
    }

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);

    int efcharge = pECharge_last[lidx];
    auto efcharged = static_cast<double>(efcharge);

    return ((efcharged * math::E_CHARGE) / pECharge_last_dt);
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getGHKI() const {
    if (pECharge_last_dt == 0) {
        return 0;
    }

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge = 0;
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        efcharge += pECharge_last[i];
    }

    auto efcharged = static_cast<double>(efcharge);

    return ((efcharged * math::E_CHARGE) / pECharge_last_dt);
}

////////////////////////////////////////////////////////////////////////////////

double Tri::computeI(double v, double dt, double simtime, double efdt) {
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (auto i: solver::ohmiccurr_local_id::range(nocs)) {
        solver::OhmicCurrdef* ocdef = patchdef()->ohmiccurrdef(i);
        // First calculate the last little bit up to the simtime
        double integral = pPoolCount[patchdef()->ohmiccurr_chanstate(i)] *
                          (simtime - pOCtime_upd[i]);
        AssertLog(integral >= 0.0);
        pOCchan_timeintg[i] += integral;
        pOCtime_upd[i] = simtime;

        // Find the mean number of channels open over the dt
        double n = pOCchan_timeintg[i] / dt;
        current += (n * ocdef->getG()) * (v - getOCerev(i));
    }
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge = 0;
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        efcharge += pECharge[i];
    }

    // The contribution from GHK charge movement.
    auto efcharged = static_cast<double>(efcharge);

    // Convert charge to coulombs and find mean current
    current += ((efcharged * math::E_CHARGE) / dt);
    resetECharge(dt, efdt, simtime);
    resetOCintegrals();

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getOhmicI(double v, double /*dt*/) const {
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (auto i: solver::ohmiccurr_local_id::range(nocs)) {
        solver::OhmicCurrdef* ocdef = patchdef()->ohmiccurrdef(i);
        // The next is ok because Patchdef returns local index
        uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(i)];
        current += (n * ocdef->getG()) * (v - getOCerev(i));
    }

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getOhmicI(solver::ohmiccurr_local_id lidx, double v, double /*dt*/) const {
    AssertLog(lidx < patchdef()->countOhmicCurrs());
    solver::OhmicCurrdef* ocdef = patchdef()->ohmiccurrdef(lidx);
    uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(lidx)];

    return (n * ocdef->getG()) * (v - getOCerev(lidx));
}

}  // namespace steps::tetexact
