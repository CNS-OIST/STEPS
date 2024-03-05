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

#include "mpi/tetvesicle/tet_rdef.hpp"

// Standard library & STL headers.
#include <algorithm>
#include <cmath>
#include <functional>

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/ghkcurr.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/raftgen.hpp"
#include "mpi/tetvesicle/raftsreac.hpp"
#include "mpi/tetvesicle/sdiff.hpp"
#include "mpi/tetvesicle/sreac.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "mpi/tetvesicle/vdepsreac.hpp"
#include "solver/ohmiccurrdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/sreacdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

TriRDEF::TriRDEF(triangle_global_id idx,
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
                 const math::point3d& position,
                 const math::point3d& trinorm,
                 int rank,
                 int host_rank)
    : pIdx(idx)
    , pPatchdef(patchdef)
    , pArea(area)
    , pInnerTet(nullptr)
    , pOuterTet(nullptr)
    , pNextTri()
    , pLengths()
    , pDist()
    , pSDiffBndDirection()
    , pECharge_last_dt(0)
    , pECharge_accum_dt(0)
    , pPosition(position)
    , pNorm(trinorm)
    , pPatchRDEF(nullptr)
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

    std::fill_n(pSDiffBndDirection, 3, false);

    pPoolOccupancy.container().resize(nspecs);
    pLastUpdate.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

TriRDEF::~TriRDEF() {
    for (auto const& i: pKProcs) {
        delete i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pPoolCount);
    util::checkpoint(cp_file, pECharge_accum);
    util::checkpoint(cp_file, pECharge_accum_dt);
    util::checkpoint(cp_file, pECharge_last);
    util::checkpoint(cp_file, pECharge_last_dt);
    util::checkpoint(cp_file, pOCtime_upd);
    util::checkpoint(cp_file, pERev);

    // NOTE not checkpointing pPoolOccupancy, pLastUpdate, pECharge, pOCchan_timeintg,
    // pAppliedRaftgens because these should be reset at time of call to checkpoint
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::restore(std::fstream& cp_file) {
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

void TriRDEF::setInnerTet(TetRDEF* t) {
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setOuterTet(TetRDEF* t) {
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setSDiffBndDirection(uint i) {
    AssertLog(i < 3);
    pSDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setNextTri(uint i, TriRDEF* t) {
    AssertLog(i <= 2);
    pNextTri[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setupKProcs(TetVesicleRDEF* tex, bool efield) {
    hasEfield = efield;

    startKProcIdx = tex->countKProcs_();

    uint j = 0;

    nKProcs = pPatchdef->countSReacs() + pPatchdef->countSurfDiffs() + patchdef()->countRaftGens() +
              patchdef()->statedef().countRaftSReacs();
    if (hasEfield) {
        nKProcs += (pPatchdef->countVDepSReacs() + pPatchdef->countGHKcurrs());
    }

    if (hostRank == myRank) {
        pKProcs.resize(nKProcs);

        // Create surface reaction kprocs
        uint nsreacs = patchdef()->countSReacs();
        for (auto i: solver::sreac_local_id::range(nsreacs)) {
            auto& srdef = patchdef()->sreacdef(i);
            auto* sr = new SReac(&srdef, this);
            AssertLog(sr != nullptr);
            pKProcs[j++] = sr;
            solver::kproc_global_id idx = tex->addKProc_(sr);
            sr->setSchedIDX(idx);
        }

        uint nsdiffs = patchdef()->countSurfDiffs();
        for (auto i: solver::surfdiff_local_id::range(nsdiffs)) {
            auto& sddef = patchdef()->surfdiffdef(i);
            auto* sd = new SDiff(&sddef, this);
            AssertLog(sd != nullptr);
            pKProcs[j++] = sd;
            solver::kproc_global_id idx = tex->addKProc_(sd);
            sd->setSchedIDX(idx);
            tex->addSDiff_(sd);
        }

        uint nrgens = patchdef()->countRaftGens();
        for (auto i: solver::raftgen_local_id::range(nrgens)) {
            auto& raftgendef = patchdef()->raftgendef(i);
            auto* raftgen = new RaftGen(&raftgendef, this);
            AssertLog(raftgen != nullptr);
            pKProcs[j++] = raftgen;
            solver::kproc_global_id idx = tex->addKProc_(raftgen);
            raftgen->setSchedIDX(idx);
        }

        uint nrsreacs = patchdef()->statedef().countRaftSReacs();
        for (auto i: solver::raftsreac_global_id::range(nrsreacs)) {
            auto& raftsreacdef = patchdef()->statedef().raftsreacdef(i);
            auto* raftsreac = new RaftSReac(&raftsreacdef, this);
            AssertLog(raftsreac != nullptr);
            pKProcs[j++] = raftsreac;
            solver::kproc_global_id idx = tex->addKProc_(raftsreac);
            raftsreac->setSchedIDX(idx);
        }

        if (efield == true) {
            uint nvdsreacs = patchdef()->countVDepSReacs();
            for (auto i: solver::vdepsreac_local_id::range(nvdsreacs)) {
                auto& vdsrdef = patchdef()->vdepsreacdef(i);
                auto* vdsr = new VDepSReac(&vdsrdef, this);
                AssertLog(vdsr != nullptr);
                pKProcs[j++] = vdsr;
                solver::kproc_global_id idx = tex->addKProc_(vdsr, true);
                vdsr->setSchedIDX(idx);
            }

            uint nghkcurrs = patchdef()->countGHKcurrs();
            for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
                auto& ghkdef = patchdef()->ghkcurrdef(i);
                auto* ghk = new GHKcurr(&ghkdef, this);
                AssertLog(ghk != nullptr);
                pKProcs[j++] = ghk;
                solver::kproc_global_id idx = tex->addKProc_(ghk, true);
                ghk->setSchedIDX(idx);
            }
        }
    } else {
        pKProcs.resize(0);
        for (uint k = 0; k < nKProcs; k++) {
            tex->addKProc_(nullptr);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setupDeps() {
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
        TriRDEF* next = nextTri(i);
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
        TetRDEF* itet = iTet();
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

        TetRDEF* otet = oTet();
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

bool TriRDEF::KProcDepSpecTet(uint kp_lidx,
                              TetRDEF* kp_container,
                              solver::spec_global_id spec_gidx) {
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

    // if kp is surface diff, which has no tet dependency, return false
    if (remain < pPatchdef->countSurfDiffs()) {
        return false;
    }
    remain -= pPatchdef->countSurfDiffs();

    // if kp is raft gen, which has no tet dependency, return false
    if (remain < pPatchdef->countRaftGens()) {
        return false;
    }
    remain -= pPatchdef->countRaftGens();

    // if kp is raft sreac
    if (remain < pPatchdef->statedef().countRaftSReacs()) {
        const auto& rsrdef = patchdef()->statedef().raftsreacdef(
            solver::raftsreac_global_id(remain));
        if (kp_container == iTet()) {
            return rsrdef.dep_I(spec_gidx) != solver::DEP_NONE;
        } else if (kp_container == oTet()) {
            return rsrdef.dep_O(spec_gidx) != solver::DEP_NONE;
        }
        return false;
    }
    remain -= pPatchdef->statedef().countRaftSReacs();

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

bool TriRDEF::KProcDepSpecTriRaftSurface(uint kp_lidx,
                                         TriRDEF* kp_container,
                                         solver::spec_global_id spec_gidx) {
    // Currently only active for RaftSReacs, because RaftDis and RaftEndocytosis
    // now have special routines and so are ignored here

    if (kp_container != this) {
        return false;
    }

    uint remain = kp_lidx;

    // if kp is surf reaction
    if (remain < pPatchdef->countSReacs()) {
        return false;
    }
    remain -= pPatchdef->countSReacs();

    // if kp is surface diff
    if (remain < pPatchdef->countSurfDiffs()) {
        return false;
    }
    remain -= pPatchdef->countSurfDiffs();

    // if kproc is raft gen
    if (remain < pPatchdef->countRaftGens()) {
        // Depends on tri species, not raft surface species
        return false;
    }
    remain -= pPatchdef->countRaftGens();

    if (remain < pPatchdef->statedef().countRaftSReacs()) {
        const auto& rsrdef = patchdef()->statedef().raftsreacdef(
            solver::raftsreac_global_id(remain));
        return rsrdef.dep_Rs(spec_gidx) != solver::DEP_NONE;
    }
    remain -= pPatchdef->statedef().countRaftSReacs();

    if (hasEfield) {
        // VDepSReac
        if (remain < pPatchdef->countVDepSReacs()) {
            return false;
        }
        remain -= pPatchdef->countVDepSReacs();

        // GHK
        if (remain < pPatchdef->countGHKcurrs()) {
            return false;
        }
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool TriRDEF::KProcDepSpecTri(uint kp_lidx,
                              TriRDEF* kp_container,
                              solver::spec_global_id spec_gidx) {
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

    // if kproc is raft gen
    if (remain < pPatchdef->countRaftGens()) {
        const auto& rgdef = patchdef()->raftgendef(solver::raftgen_local_id(remain));
        return rgdef.dep_S(spec_gidx) != solver::DEP_NONE;
    }
    remain -= pPatchdef->countRaftGens();

    if (remain < pPatchdef->statedef().countRaftSReacs()) {
        const auto& rsrdef = patchdef()->statedef().raftsreacdef(
            solver::raftsreac_global_id(remain));
        return rsrdef.dep_S(spec_gidx) != solver::DEP_NONE;
    }
    remain -= pPatchdef->statedef().countRaftSReacs();

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

void TriRDEF::reset() {
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

    clearRaftProxyrefs();
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::resetECharge(double dt, double efdt) {
    const uint nghkcurrs = pPatchdef->countGHKcurrs();
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        pECharge_accum[i] += pECharge[i];
    }
    pECharge_accum_dt += dt;

    if (pECharge_accum_dt >= efdt) {
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

void TriRDEF::resetOCintegrals() {
    std::fill(pOCchan_timeintg.begin(), pOCchan_timeintg.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::incECharge(solver::ghkcurr_local_id lidx, int charge) {
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);
    pECharge[lidx] += charge;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setCount(solver::spec_local_id lidx, uint count, double period) {
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

void TriRDEF::incCount(solver::spec_local_id lidx, int inc, double period, bool local_change) {
    // can change count locally or remotely, if change locally, check if the
    // change require sync, if change remotely, register the change in pRDEF, sync
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
        bufferLocations[lidx] = pRDEF->registerRemoteMoleculeChange_(
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

void TriRDEF::setOCchange(solver::ohmiccurr_local_id oclidx,
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

void TriRDEF::setOCerev(solver::ohmiccurr_local_id oclidx, double erev) {
    AssertLog(oclidx < patchdef()->countOhmicCurrs());
    pERev[oclidx] = erev;
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getOCerev(solver::ohmiccurr_local_id oclidx) const {
    auto it = pERev.find(oclidx);
    if (it != pERev.end()) {
        return it->second;
    } else {
        return patchdef()->ohmiccurrdef(oclidx).getERev();
    }
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setClamped(solver::spec_local_id lidx, bool clamp) {
    if (clamp == true) {
        pPoolFlags[lidx] |= CLAMPED;
    } else {
        pPoolFlags[lidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

int TriRDEF::getTriDirection(triangle_global_id tidx) {
    for (uint i = 0; i < 3; i++) {
        if (pTris[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

SReac& TriRDEF::sreac(solver::sreac_local_id lidx) const {
    AssertLog(lidx < patchdef()->countSReacs());
    return *dynamic_cast<SReac*>(pKProcs[lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

SDiff& TriRDEF::sdiff(solver::surfdiff_local_id lidx) const {
    AssertLog(lidx < patchdef()->countSurfDiffs());
    return *dynamic_cast<SDiff*>(pKProcs[patchdef()->countSReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

RaftGen& TriRDEF::raftgen(solver::raftgen_local_id lidx) const {
    AssertLog(lidx < patchdef()->countRaftGens());
    return *dynamic_cast<RaftGen*>(
        pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

VDepSReac& TriRDEF::vdepsreac(solver::vdepsreac_local_id lidx) const {
    AssertLog(lidx < patchdef()->countVDepSReacs());
    return *dynamic_cast<VDepSReac*>(
        pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() +
                patchdef()->countRaftGens() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

GHKcurr& TriRDEF::ghkcurr(solver::ghkcurr_local_id lidx) const {
    AssertLog(lidx < patchdef()->countGHKcurrs());
    return *dynamic_cast<GHKcurr*>(
        pKProcs[patchdef()->countSReacs() + patchdef()->countSurfDiffs() +
                patchdef()->countRaftGens() + patchdef()->countVDepSReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getGHKI(solver::ghkcurr_local_id lidx, double dt) const {
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    AssertLog(lidx < nghkcurrs);

    int efcharge = pECharge_last[lidx];
    auto efcharged = static_cast<double>(efcharge);

    return (efcharged * math::E_CHARGE) / dt;
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getGHKI(double dt) const {
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge = 0;
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        efcharge += pECharge_last[i];
    }

    auto efcharged = static_cast<double>(efcharge);

    return (efcharged * math::E_CHARGE) / dt;
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::computeI(double v, double dt, double simtime, double efdt) {
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
    uint nghkcurrs = pPatchdef->countGHKcurrs();
    int efcharge = 0;
    for (auto i: solver::ghkcurr_local_id::range(nghkcurrs)) {
        efcharge += pECharge[i];
    }

    // The contribution from GHK charge movement.
    auto efcharged = static_cast<double>(efcharge);

    // Convert charge to coulombs and find mean current
    current += ((efcharged * math::E_CHARGE) / dt);
    resetECharge(dt, efdt);
    resetOCintegrals();

    return current;
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getOhmicI(double v, double /*dt*/) const {
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

double TriRDEF::getOhmicI(solver::ohmiccurr_local_id lidx, double v, double /*dt*/) const {
    AssertLog(lidx < patchdef()->countOhmicCurrs());
    const auto& ocdef = patchdef()->ohmiccurrdef(lidx);
    uint n = pPoolCount[patchdef()->ohmiccurr_chanstate(lidx)];
    return (n * ocdef.getG()) * (v - getOCerev(lidx));
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::createRaftProxyref(solver::Raftdef* raftdef,
                                 solver::raft_individual_id raft_unique_id) {
    if (pRaftProxyrefs.find(raft_unique_id) != pRaftProxyrefs.end()) {
        ProgErrLog("RaftProxy already assigned to Tri.\n");
    }

    auto* raftproxy = new RaftProxy(raftdef, this, raft_unique_id);

    pRaftProxyrefs[raft_unique_id] = raftproxy;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::clearRaftProxyrefs() {
    for (auto& raft_proxy: pRaftProxyrefs) {
        delete raft_proxy.second;
    }
    pRaftProxyrefs.clear();
}

////////////////////////////////////////////////////////////////////////////////

RaftProxy* TriRDEF::getRaftProxyref(solver::raft_individual_id raft_uidx) {
    if (pRaftProxyrefs.find(raft_uidx) == pRaftProxyrefs.end()) {
        ProgErrLog("RaftProxy unknown in Tri.\n");
    }

    return pRaftProxyrefs[raft_uidx];
}

////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool TriRDEF::getInHost() const {
    return hostRank == myRank;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setHost(int host, int rank) {
    hostRank = host;
    myRank = rank;
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setSolverRDEF(TetVesicleRDEF* solver) {
    pRDEF = solver;
}

////////////////////////////////////////////////////////////////////////////////

TetVesicleRDEF* TriRDEF::solverRDEF() const {
    return pRDEF;
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getPoolOccupancy(solver::spec_local_id lidx) {
    return pPoolOccupancy.at(lidx);
}

////////////////////////////////////////////////////////////////////////////////

double TriRDEF::getLastUpdate(solver::spec_local_id lidx) {
    return pLastUpdate.at(lidx);
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::resetPoolOccupancy() {
    std::fill(pPoolOccupancy.begin(), pPoolOccupancy.end(), 0.0);
    std::fill(pLastUpdate.begin(), pLastUpdate.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::setupBufferLocations() {
    uint nspecs = pPatchdef->countSpecs();
    bufferLocations.container().assign(nspecs, std::numeric_limits<uint>::max());
}

////////////////////////////////////////////////////////////////////////////////

void TriRDEF::applyRaftGen(solver::RaftGendef* rgdef, double period) {
    // Ad the raftgen to the list to apply at next update,
    // and remove the species from tri pools.

    pAppliedRaftgens[rgdef->gidx()] += 1;

    // This uses global indices

    const auto& lhs_s_vec = rgdef->lhs_S();

    for (auto sg: solver::spec_global_id::range(lhs_s_vec.size())) {
        uint lhs = lhs_s_vec[sg];
        if (lhs == 0) {
            continue;
        }
        // Set the raft count- now needs to be done during patch application
        // patch->setSingleRaftSpecCount(rgidx, raft_unique_index, sg, lhs);

        // Set the tri count
        solver::spec_local_id spec_lidx = patchdef()->specG2L(sg);

        AssertLog(spec_lidx.valid());

        if (clamped(spec_lidx) == true) {
            continue;
        }
        uint prev_count = pools()[spec_lidx];

        int new_count = static_cast<int>(prev_count) - static_cast<int>(lhs);

        AssertLog(new_count >= 0);

        setCount(spec_lidx, static_cast<uint>(new_count), period);
    }
}

}  // namespace steps::mpi::tetvesicle
