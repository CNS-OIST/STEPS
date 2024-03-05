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

#include <algorithm>
#include <cmath>
#include <functional>
#include <iostream>


#include "ghkcurr.hpp"
#include "sdiff.hpp"
#include "sreac.hpp"
#include "tet.hpp"
#include "tetopsplit.hpp"
#include "vdepsreac.hpp"

#include "math/constants.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/sreacdef.hpp"

#include "util/checkpointing.hpp"
#include "util/error.hpp"

namespace steps::mpi::tetopsplit {


////////////////////////////////////////////////////////////////////////////////

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
         triangle_global_id tri2,
         int rank,
         int host_rank)
    : pIdx(idx)
    , pPatchdef(patchdef)
    , pArea(area)
    , pLengths()
    , pDist()
    , pNextTri()
    , pECharge_last_dt(0)
    , pECharge_accum_dt(0)
    , hostRank(host_rank)
    , myRank(rank) {
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
    pPoolCount.container().resize(nspecs);
    pPoolFlags.container().resize(nspecs);

    uint nghkcurrs = pPatchdef->countGHKcurrs();
    pECharge.container().resize(nghkcurrs);
    pECharge_last.container().resize(nghkcurrs);
    pECharge_accum.container().resize(nghkcurrs);

    uint nohmcurrs = pPatchdef->countOhmicCurrs();
    pOCchan_timeintg.container().resize(nohmcurrs);
    pOCtime_upd.container().resize(nohmcurrs);

    pPoolOccupancy.container().resize(nspecs);
    pLastUpdate.container().resize(nspecs);

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
    // NOTE not checkpointing pPoolOccupancy, pLastUpdate, pECharge, pOCchan_timeintg because these
    // should be reset at time of call to checkpoint
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
    AssertLog(i <= 2);

    pNextTri[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupKProcs(TetOpSplitP* tex, bool efield) {
    hasEfield = efield;
    startKProcIdx = tex->countKProcs();
    uint j = 0;

    nKProcs = pPatchdef->countSReacs() + pPatchdef->countSurfDiffs();
    if (hasEfield) {
        // pPatchdef->countGHKcurrs());
        nKProcs += (pPatchdef->countVDepSReacs() + pPatchdef->countGHKcurrs());
    }

    if (hostRank == myRank) {
        pKProcs.resize(nKProcs);

        // Create surface reaction kprocs
        uint nsreacs = patchdef()->countSReacs();
        for (auto i: solver::sreac_local_id::range(nsreacs)) {
            auto& srdef = patchdef()->sreacdef(i);
            auto* sr = new SReac(&srdef, this);
            pKProcs[j++] = sr;
            solver::kproc_global_id idx = tex->addKProc(sr);
            sr->setSchedIDX(idx);
        }

        uint nsdiffs = patchdef()->countSurfDiffs();
        for (auto i: solver::surfdiff_local_id::range(nsdiffs)) {
            auto& sddef = patchdef()->surfdiffdef(i);
            auto* sd = new SDiff(&sddef, this);
            pKProcs[j++] = sd;
            solver::kproc_global_id idx = tex->addKProc(sd);
            sd->setSchedIDX(idx);
            tex->addSDiff(sd);
        }

        // REMINDER: hasEfield not ready
        if (hasEfield) {
            uint nvdsreacs = patchdef()->countVDepSReacs();
            for (auto i: solver::vdepsreac_local_id::range(nvdsreacs)) {
                auto& vdsrdef = patchdef()->vdepsreacdef(i);
                auto* vdsr = new VDepSReac(&vdsrdef, this);
                pKProcs[j++] = vdsr;
                solver::kproc_global_id idx = tex->addKProc(vdsr, true);
                vdsr->setSchedIDX(idx);
            }

            uint nghkcurrs = patchdef()->countGHKcurrs();
            for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
                auto& ghkdef = patchdef()->ghkcurrdef(i);
                auto* ghk = new GHKcurr(&ghkdef, this);
                pKProcs[j++] = ghk;
                solver::kproc_global_id idx = tex->addKProc(ghk, true);
                ghk->setSchedIDX(idx);
            }
        }
    } else {
        pKProcs.resize(0);
        for (uint k = 0; k < nKProcs; k++) {
            tex->addKProc(nullptr);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupDeps() {
    if (myRank != hostRank) {
        return;
    }
    for (auto& kp: pKProcs) {
        kp->setupDeps();
    }

    bool has_remote_neighbors = false;
    uint nspecs = patchdef()->countSpecs();
    for (uint i = 0; i < 3; ++i) {
        // Fetch next triangles, if it exists.
        Tri* next = nextTri(i);
        if (next == nullptr) {
            continue;
        }
        if (next->getHost() != getHost()) {
            has_remote_neighbors = true;
            break;
        }
    }
    if (!has_remote_neighbors) {
        localSpecUpdKProcs.container().clear();
        return;
    }
    localSpecUpdKProcs.container().resize(nspecs);

    for (auto slidx: solver::spec_local_id::range(nspecs)) {
        solver::spec_global_id sgidx = patchdef()->specL2G(slidx);
        // check dependency of kprocs in the same tri
        uint nkprocs = countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            // Check locally.
            if (KProcDepSpecTri(sk, this, sgidx)) {
                localSpecUpdKProcs[slidx].push_back(getKProc(sk));
            }
        }

        // Check the neighbouring tetrahedrons.
        WmVol* itet = iTet();
        if (itet != nullptr) {
            if (getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << idx() << " and its compartment tetrahedron "
                   << itet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = itet->countKProcs();
            for (uint k = 0; k < nkprocs; k++) {
                if (itet->KProcDepSpecTri(k, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(itet->getKProc(k));
                }
            }
        }

        WmVol* otet = oTet();
        if (otet != nullptr) {
            if (getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << idx() << " and its compartment tetrahedron "
                   << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = otet->countKProcs();
            for (uint k = 0; k < nkprocs; k++) {
                if (otet->KProcDepSpecTri(k, this, sgidx)) {
                    localSpecUpdKProcs[slidx].push_back(otet->getKProc(k));
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Tri::KProcDepSpecTet(uint kp_lidx, WmVol* kp_container, solver::spec_global_id spec_gidx) {
    // if kp is surf reaction
    uint remain = kp_lidx;
    if (remain < pPatchdef->countSReacs()) {
        const auto& srdef = patchdef()->sreacdef(solver::sreac_local_id(remain));
        if (kp_container == iTet()) {
            return srdef.dep_I(spec_gidx) != solver::DEP_NONE;
        } else if (kp_container == oTet()) {
            return srdef.dep_O(spec_gidx) != solver::DEP_NONE;
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
    if (hasEfield) {
        // VDepSReac
        if (remain < pPatchdef->countVDepSReacs()) {
            const auto& vdsrdef = patchdef()->vdepsreacdef(solver::vdepsreac_local_id(remain));
            if (kp_container == iTet()) {
                return vdsrdef.dep_I(spec_gidx) != solver::DEP_NONE;
            } else if (kp_container == oTet()) {
                return vdsrdef.dep_O(spec_gidx) != solver::DEP_NONE;
            }
            return false;
        }
        remain -= pPatchdef->countVDepSReacs();

        // GHK
        if (remain < pPatchdef->countGHKcurrs()) {
            const auto& ghkdef = patchdef()->ghkcurrdef(solver::ghkcurr_local_id(remain));
            if (kp_container == iTet()) {
                return ghkdef.dep_v(spec_gidx) != solver::DEP_NONE;
            } else if (kp_container == oTet()) {
                if (ghkdef.voconc() < 0.0) {
                    return ghkdef.dep_v(spec_gidx) != solver::DEP_NONE;
                } else {
                    return false;
                }
            }

            return false;
        }
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool Tri::KProcDepSpecTri(uint kp_lidx, Tri* kp_container, solver::spec_global_id spec_gidx) {
    if (kp_container != this) {
        return false;
    }

    // if kp is surf reaction
    uint remain = kp_lidx;
    if (remain < pPatchdef->countSReacs()) {
        const auto& srdef = patchdef()->sreacdef(solver::sreac_local_id(remain));
        return srdef.dep_S(spec_gidx) != solver::DEP_NONE;
    }
    remain -= pPatchdef->countSReacs();

    // if kp is surface diff
    if (remain < pPatchdef->countSurfDiffs()) {
        const auto& sddef = patchdef()->surfdiffdef(solver::surfdiff_local_id(remain));
        return spec_gidx == sddef.lig();
    }

    remain -= pPatchdef->countSurfDiffs();

    // REMINDER: hasEfield not ready
    if (hasEfield) {
        // VDepSReac
        if (remain < pPatchdef->countVDepSReacs()) {
            const auto& vdsrdef = patchdef()->vdepsreacdef(solver::vdepsreac_local_id(remain));
            return vdsrdef.dep_S(spec_gidx) != solver::DEP_NONE;
        }
        remain -= pPatchdef->countVDepSReacs();

        // GHK
        if (remain < pPatchdef->countGHKcurrs()) {
            const auto& ghkdef = patchdef()->ghkcurrdef(solver::ghkcurr_local_id(remain));
            return ghkdef.dep(spec_gidx) != solver::DEP_NONE;
        }
    }

    AssertLog(false);
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
    for (auto i: pECharge.range()) {
        pECharge_accum[i] += pECharge[i];
    }
    pECharge_accum_dt += dt;

    if (pECharge_accum_dt >= efdt or
        (efdt - pECharge_accum_dt) <= std::numeric_limits<double>::epsilon() * t * 8) {
        // Swap arrays
        std::swap(pECharge_last, pECharge_accum);

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
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);
    pECharge[lidx] += charge;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setCount(solver::spec_local_id lidx, uint count, double period) {
    // tri count doesn't care about global
    AssertLog(lidx < patchdef()->countSpecs());
    double oldcount = pPoolCount[lidx];
    pPoolCount[lidx] = count;
    if (period == 0.0) {
        return;
    }
    // Count has changed,
    double lastupdate = pLastUpdate[lidx];
    AssertLog(period >= lastupdate);
    pPoolOccupancy[lidx] += oldcount * (period - lastupdate);

    pLastUpdate[lidx] = period;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::incCount(solver::spec_local_id lidx, int inc, double period, bool local_change) {
    // can change count locally or remotely, if change locally, check if the
    // change require sync, if change remotely, register the change in pSol, sync
    // of remote change will be register when the change is applied in remote host
    // via this function
    AssertLog(lidx < patchdef()->countSpecs());

    // count changed by diffusion
    if (hostRank != myRank && !local_change) {
        if (inc <= 0) {
            std::ostringstream os;
            os << "Try to change molecule " << lidx << " by " << inc << "\n";
            os << "Fail because molecule change of receiving end should always be "
                  "non-negative.\n";
            ProgErrLog(os.str());
        }
        bufferLocations[lidx] = pSol->registerRemoteMoleculeChange(
            hostRank, bufferLocations[lidx], SUB_TRI, pIdx.get(), lidx, inc);
    }
    // local change by reac or diff
    else {
        // need to check if the change costs negative molecule count
        double oldcount = pPoolCount[lidx];
        AssertLog(oldcount + inc >= 0.0);
        pPoolCount[lidx] += inc;

        if (period == 0.0 || local_change) {
            return;
        }
        // Count has changed,
        double lastupdate = pLastUpdate[lidx];
        AssertLog(period >= lastupdate);
        pPoolOccupancy[lidx] += oldcount * (period - lastupdate);

        pLastUpdate[lidx] = period;
    }
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
        return patchdef()->ohmiccurrdef(oclidx).getERev();
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
int Tri::getTriDirection(triangle_global_id tidx) const {
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

    return (efcharged * math::E_CHARGE) / pECharge_last_dt;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getGHKI() const {
    if (pECharge_last_dt == 0) {
        return 0;
    }

    int efcharge = std::accumulate(pECharge_last.begin(), pECharge_last.end(), 0);

    auto efcharged = static_cast<double>(efcharge);

    return (efcharged * math::E_CHARGE) / pECharge_last_dt;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::computeI(double v, double dt, double simtime, double efdt) {
    double current = 0.0;
    uint nocs = patchdef()->countOhmicCurrs();
    for (auto i: solver::ohmiccurr_local_id::range(nocs)) {
        const auto& ocdef = patchdef()->ohmiccurrdef(i);
        // First calculate the last little bit up to the simtime
        double integral = pPoolCount[patchdef()->ohmiccurr_chanstate(i)] *
                          (simtime - pOCtime_upd[i]);
        AssertLog(integral >= 0.0);
        pOCchan_timeintg[i] += integral;
        pOCtime_upd[i] = simtime;

        // Find the mean number of channels open over the dt
        double n = pOCchan_timeintg[i] / dt;
        current += (n * ocdef.getG()) * (v - getOCerev(i));
    }

    int efcharge = std::accumulate(pECharge.begin(), pECharge.end(), 0);

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
        const auto& ocdef = patchdef()->ohmiccurrdef(i);
        // The next is ok because Patchdef returns local index
        uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(i)];
        current += (n * ocdef.getG()) * (v - getOCerev(i));
    }

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getOhmicI(solver::ohmiccurr_local_id lidx, double v, double /*dt*/) const {
    AssertLog(lidx < patchdef()->countOhmicCurrs());
    const auto& ocdef = patchdef()->ohmiccurrdef(lidx);
    uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(lidx)];

    return (n * ocdef.getG()) * (v - getOCerev(lidx));
}

////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool Tri::getInHost() const {
    return hostRank == myRank;
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setHost(int host, int rank) {
    hostRank = host;
    myRank = rank;
}
////////////////////////////////////////////////////////////////////////////////

void Tri::setSolver(mpi::tetopsplit::TetOpSplitP* sol) {
    pSol = sol;
}

////////////////////////////////////////////////////////////////////////////////

TetOpSplitP* Tri::solver() {
    return pSol;
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getPoolOccupancy(solver::spec_local_id lidx) {
    AssertLog(lidx < patchdef()->countSpecs());

    return pPoolOccupancy[lidx];
}

////////////////////////////////////////////////////////////////////////////////

double Tri::getLastUpdate(solver::spec_local_id lidx) {
    AssertLog(lidx < patchdef()->countSpecs());

    return pLastUpdate[lidx];
}

////////////////////////////////////////////////////////////////////////////////

void Tri::resetPoolOccupancy() {
    std::fill(pPoolOccupancy.begin(), pPoolOccupancy.end(), 0.0);
    std::fill(pLastUpdate.begin(), pLastUpdate.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void Tri::repartition(TetOpSplitP* tex, int rank, int host_rank) {
    myRank = rank;
    hostRank = host_rank;

    // Delete reaction rules.
    auto e = pKProcs.end();
    for (auto i = pKProcs.begin(); i != e; ++i) {
        delete *i;
    }

    setupKProcs(tex);

    localSpecUpdKProcs.container().clear();
    bufferLocations.container().clear();
}

////////////////////////////////////////////////////////////////////////////////

void Tri::setupBufferLocations() {
    uint nspecs = pPatchdef->countSpecs();
    bufferLocations.container().assign(nspecs, std::numeric_limits<uint>::max());
}

}  // namespace steps::mpi::tetopsplit
