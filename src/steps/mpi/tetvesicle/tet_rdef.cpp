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
#include <functional>

// STEPS headers.
#include "mpi/tetvesicle/diff.hpp"
#include "mpi/tetvesicle/exocytosis.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/reac.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "mpi/tetvesicle/vesbind.hpp"
#include "mpi/tetvesicle/vesreac.hpp"
#include "solver/compdef.hpp"
#include "solver/diffdef.hpp"
#include "solver/reacdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

TetRDEF::TetRDEF(tetrahedron_global_id idx,
                 solver::Compdef* cdef,
                 double vol,
                 double a0,
                 double a1,
                 double a2,
                 double a3,
                 double d0,
                 double d1,
                 double d2,
                 double d3,
                 tetrahedron_global_id tet0,
                 tetrahedron_global_id tet1,
                 tetrahedron_global_id tet2,
                 tetrahedron_global_id tet3,
                 const math::point3d& /*baryc*/,
                 int rank,
                 int host_rank)
    : pIdx(idx)
    , pCompdef(cdef)
    , pVol(vol)
    , pOverlap(0.0)
    , pCompRDEF(nullptr)
    , pNextTet()
    , pAreas()
    , pDist()
    , myRank(rank)
    , hostRank(host_rank) {
    AssertLog(pCompdef != nullptr);
    AssertLog(pVol > 0.0);

    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount.container().resize(nspecs);
    pPoolFlags.container().resize(nspecs);
    pPoolOccupancy.container().resize(nspecs);
    pLastUpdate.container().resize(nspecs);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i = 0; i <= 3; ++i) {
        pNextTet[i] = nullptr;
        pNextTris[i] = nullptr;
    }

    pTets[0] = tet0;
    pTets[1] = tet1;
    pTets[2] = tet2;
    pTets[3] = tet3;

    pAreas[0] = a0;
    pAreas[1] = a1;
    pAreas[2] = a2;
    pAreas[3] = a3;

    pDist[0] = d0;
    pDist[1] = d1;
    pDist[2] = d2;
    pDist[3] = d3;

    pDiffBndDirection.fill(false);

    nKProcs = compdef()->countDiffs() + compdef()->countReacs() + compdef()->countVesBinds() +
              compdef()->statedef().countVesSReacs() + compdef()->statedef().countExocytosis();
}

////////////////////////////////////////////////////////////////////////////////

TetRDEF::~TetRDEF() {
    // Delete kprocs.
    for (auto const& i: pKProcs) {
        delete i;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolFlags);
    util::checkpoint(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolFlags);
    util::restore(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setNextTet(uint i, TetRDEF* t) {
    // Now adding all tets, even those from other compartments, due to the
    // diffusion boundaries
    pNextTet[i] = t;

    pNextTris[i] = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setDiffBndDirection(uint i) {
    AssertLog(i < 4);

    pDiffBndDirection[i] = true;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setNextTri(uint i, TriRDEF* t) {
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    pNextTet[i] = nullptr;
    pNextTris[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::reset() {
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
    std::fill(pPoolFlags.begin(), pPoolFlags.end(), 0);

    for (auto const& kproc: pKProcs) {
        kproc->reset();
    }

    pOverlap = 0.0;

    clearVesProxyrefs();
}

////////////////////////////////////////////////////////////////////////////////
// Add a vesicle reference (vesicle overlaps this tet)
void TetRDEF::createVesProxyref(solver::Vesicledef* vesdef,
                                solver::vesicle_individual_id ves_unique_id,
                                math::position_abs ves_pos,
                                bool contains_link) {
    if (pVesProxyrefs.find(ves_unique_id) != pVesProxyrefs.end()) {
        ProgErrLog("VesProxy already assigned to Tet.\n");
    }

    auto* vesproxy = new VesProxy(vesdef, this, ves_unique_id, ves_pos, contains_link);

    pVesProxyrefs[ves_unique_id] = vesproxy;
}

////////////////////////////////////////////////////////////////////////////////

VesProxy* TetRDEF::getVesProxyref(solver::vesicle_individual_id ves_uidx) {
    if (pVesProxyrefs.find(ves_uidx) == pVesProxyrefs.end()) {
        ProgErrLog("VesProxy unknown in Tet.\n");
    }

    return pVesProxyrefs[ves_uidx];
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::clearVesProxyrefs() {
    for (auto& ves_proxy: pVesProxyrefs) {
        delete ves_proxy.second;
    }
    pVesProxyrefs.clear();
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setupKProcs(TetVesicleRDEF* tex) {
    startKProcIdx = tex->countKProcs_();

    uint j = 0;

    uint nreacs = compdef()->countReacs();
    uint ndiffs = compdef()->countDiffs();
    uint vesbinds = compdef()->countVesBinds();

    // Used to be vesicle-based kprocs, now tet-based for MPISTEPS
    uint nvsreacs = compdef()->statedef().countVesSReacs();
    uint nexos = compdef()->statedef().countExocytosis();

    AssertLog(nKProcs == nreacs + ndiffs + vesbinds + nvsreacs + nexos);

    // if in host create KProc
    if (hostRank == myRank) {
        pKProcs.resize(nKProcs);

        for (auto i: solver::reac_local_id::range(nreacs)) {
            auto& rdef = compdef()->reacdef(i);
            Reac* r = new Reac(&rdef, this);
            kprocs()[j++] = r;
            solver::kproc_global_id idx = tex->addKProc_(r);
            r->setSchedIDX(idx);
        }

        for (auto i: solver::diff_local_id::range(ndiffs)) {
            auto& ddef = compdef()->diffdef(i);
            Diff* d = new Diff(&ddef, this);
            kprocs()[j++] = d;
            solver::kproc_global_id idx = tex->addKProc_(d);
            d->setSchedIDX(idx);
            tex->addDiff_(d);
        }

        for (auto i: solver::vesbind_local_id::range(vesbinds)) {
            auto& vbdef = compdef()->vesbinddef(i);
            auto* vb = new VesBind(&vbdef, this);
            kprocs()[j++] = vb;
            solver::kproc_global_id idx = tex->addKProc_(vb);
            vb->setSchedIDX(idx);
        }

        // Used to be vesicle-based kprocs, now tet-based for MPISTEPS
        for (auto i: solver::vessreac_global_id::range(nvsreacs)) {
            auto& vsrdef = compdef()->statedef().vessreacdef(i);
            auto* vr = new VesReac(&vsrdef, this);
            kprocs()[j++] = vr;
            solver::kproc_global_id idx = tex->addKProc_(vr);
            vr->setSchedIDX(idx);
        }
        for (auto i: solver::exocytosis_global_id::range(nexos)) {
            auto& exodef = compdef()->statedef().exocytosisdef(i);
            auto* exo = new Exocytosis(&exodef, this);
            kprocs()[j++] = exo;
            solver::kproc_global_id idx = tex->addKProc_(exo);
            exo->setSchedIDX(idx);
        }

    }
    // else just record the idx. TODO is this needed??
    else {
        kprocs().resize(0);

        for (uint i = 0; i < nKProcs; ++i) {
            tex->addKProc_(nullptr);
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setupDeps() {
    if (myRank != hostRank) {
        return;
    }

    for (auto& kp: pKProcs) {
        kp->setupDeps();
    }

    bool has_remote_neighbors = false;
    const uint nspecs = compdef()->countSpecs();
    for (uint i = 0; i < 4; ++i) {
        // Fetch next tetrahedron, if it exists.
        TetRDEF* next = nextTet(i);
        if (next == nullptr) {
            continue;
        }
        if (next->getHost() != getHost()) {
            has_remote_neighbors = true;
            break;
        }
    }

    if (has_remote_neighbors == false) {
        localSpecUpdKProcs.clear();
        return;
    }

    localSpecUpdKProcs.resize(nspecs);
    for (auto slidx: solver::spec_local_id::range(nspecs)) {
        solver::spec_global_id sgidx = compdef()->specL2G(slidx);
        // search dependency for kprocs in this tet
        uint nkprocs = countKProcs();

        for (uint k = 0; k < nkprocs; k++) {
            if (KProcDepSpecTet(k, this, sgidx)) {
                localSpecUpdKProcs[slidx.get()].push_back(getKProc(k));
            }
        }

        // search dependency for kprocs in neighboring tris
        for (uint i = 0; i < 4; ++i) {
            TriRDEF* next = nextTri(i);
            if (next == nullptr) {
                continue;
            }

            // next tri has to be in the same host to prevent
            // cross process surface reaction
            if (next->getHost() != getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next->idx() << " and its compartment tetrahedron "
                   << idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            nkprocs = next->countKProcs();
            for (uint sk = 0; sk < nkprocs; sk++) {
                if (next->KProcDepSpecTet(sk, this, sgidx)) {
                    localSpecUpdKProcs[slidx.get()].push_back(next->getKProc(sk));
                }
            }
        }
    }
}

////////////////////////////////////////////////////////////////////////////////

bool TetRDEF::KProcDepSpecTetVesSurface(uint kp_lidx,
                                        TetRDEF* kp_container,
                                        solver::spec_global_id spec_gidx) {
    if (kp_container != this) {
        return false;
    }
    uint remain = kp_lidx;

    // if kp is reaction
    if (remain < compdef()->countReacs()) {
        return false;
    }
    remain -= compdef()->countReacs();

    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        return false;
    }
    remain -= compdef()->countDiffs();

    // if kp is vesbind
    if (remain < compdef()->countVesBinds()) {
        const auto& vbdef = compdef()->vesbinddef(solver::vesbind_local_id(remain));
        if (spec_gidx == vbdef.getSpec1gidx() || spec_gidx == vbdef.getSpec2gidx()) {
            return true;
        }
        if ((vbdef.vdep1(spec_gidx) != 0) || (vbdef.vdep2(spec_gidx) != 0)) {
            return true;
        }
        return false;
    }
    remain -= compdef()->countVesBinds();

    // if kp is vessreac
    if (remain < compdef()->statedef().countVesSReacs()) {
        const auto& vsrdef = compdef()->statedef().vessreacdef(solver::vessreac_global_id(remain));
        return vsrdef.dep_V(spec_gidx) != solver::DEP_NONE;
    }
    remain -= compdef()->statedef().countVesSReacs();

    // if kp is exocytosis
    if (remain < compdef()->statedef().countExocytosis()) {
        const auto& exodef = compdef()->statedef().exocytosisdef(
            solver::exocytosis_global_id(remain));
        return exodef.dep_V(spec_gidx) != solver::DEP_NONE;
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool TetRDEF::KProcDepSpecTet(uint kp_lidx,
                              TetRDEF* kp_container,
                              solver::spec_global_id spec_gidx) {
    if (kp_container != this) {
        return false;
    }

    // if kp is reaction
    uint remain = kp_lidx;
    if (remain < compdef()->countReacs()) {
        const auto& rdef = compdef()->reacdef(solver::reac_local_id(remain));
        return rdef.dep(spec_gidx) != solver::DEP_NONE;
    }
    remain -= compdef()->countReacs();

    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        return spec_gidx == compdef()->diffdef(solver::diff_local_id(remain)).lig();
    }
    remain -= compdef()->countDiffs();

    // if kp is vesbind
    if (remain < compdef()->countVesBinds()) {
        // Doesn't depend on tet species persay, depends on vesicle surface
        return false;
    }
    remain -= compdef()->countVesBinds();

    // if kp is vessreac
    if (remain < compdef()->statedef().countVesSReacs()) {
        const auto& vsrdef = compdef()->statedef().vessreacdef(solver::vessreac_global_id(remain));
        return vsrdef.dep_O(spec_gidx) != solver::DEP_NONE;
    }
    remain -= compdef()->statedef().countVesSReacs();

    // if kp is exocytosis
    if (remain < compdef()->statedef().countExocytosis()) {
        return false;  // no tet dependency for exocytosis- it is vessurface based
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

bool TetRDEF::KProcDepLinkSpecTetVesSurface(uint kp_lidx,
                                            TetRDEF* kp_container,
                                            solver::linkspec_global_id linkspec_gidx) {
    if (kp_container != this) {
        return false;
    }

    // if kp is reaction
    uint remain = kp_lidx;
    if (remain < compdef()->countReacs()) {
        return false;  // no linkspec dependency for reac
    }
    remain -= compdef()->countReacs();

    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        return false;  // no linkspec dependency for diff
    }
    remain -= compdef()->countDiffs();

    // if kp is vesbind
    if (remain < compdef()->countVesBinds()) {
        const auto& vbdef = compdef()->vesbinddef(solver::vesbind_local_id(remain));
        if ((vbdef.ldep1(linkspec_gidx) != 0) || (vbdef.ldep2(linkspec_gidx) != 0)) {
            return true;
        }
        return false;
    }
    remain -= compdef()->countVesBinds();

    // if kp is vessreac
    if (remain < compdef()->statedef().countVesSReacs()) {
        const auto& vsrdef = compdef()->statedef().vessreacdef(solver::vessreac_global_id(remain));
        return vsrdef.dep_L(linkspec_gidx) != solver::DEP_NONE;
    }
    remain -= compdef()->statedef().countVesSReacs();

    // if kp is exocytosis
    if (remain < compdef()->statedef().countExocytosis()) {
        return false;
    }

    AssertLog(false);
}
////////////////////////////////////////////////////////////////////////////////

bool TetRDEF::KProcDepSpecTri(uint kp_lidx,
                              TriRDEF* kp_container,
                              solver::spec_global_id spec_gidx) const {
    if (std::find(nexttriBegin(), nexttriEnd(), kp_container) == nexttriEnd()) {
        return false;
    }

    // if kp is reaction
    uint remain = kp_lidx;
    if (remain < compdef()->countReacs()) {
        return false;
    }
    remain -= compdef()->countReacs();

    // if kp is  diff
    if (remain < compdef()->countDiffs()) {
        return false;
    }
    remain -= compdef()->countDiffs();

    // if kp is vesbind
    if (remain < compdef()->countVesBinds()) {
        return false;
    }
    remain -= compdef()->countVesBinds();

    // if kp is vessreac
    if (remain < compdef()->statedef().countVesSReacs()) {
        const auto& vsrdef = compdef()->statedef().vessreacdef(solver::vessreac_global_id(remain));
        return vsrdef.dep_S(spec_gidx) != solver::DEP_NONE;
    }
    remain -= compdef()->statedef().countVesSReacs();

    // if kp is exocytosis
    if (remain < compdef()->statedef().countExocytosis()) {
        return false;
    }

    AssertLog(false);
}

////////////////////////////////////////////////////////////////////////////////

double TetRDEF::conc(solver::spec_global_id gidx) const {
    solver::spec_local_id lspidx = compdef()->specG2L(gidx);
    double n = pPoolCount[lspidx];
    return n / (1.0e3 * pVol * math::AVOGADRO);
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setCount(solver::spec_local_id lidx, uint count, double period) {
    // Count has changed, need to correct pool factor

    // This function is used for updates that do not require remote sync,
    // such as user setCount or reaction

    AssertLog(lidx < compdef()->countSpecs());
    uint oldcount = pPoolCount[lidx];
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

void TetRDEF::incCount(solver::spec_local_id lidx, int inc, double period, bool local_change) {
    AssertLog(lidx < compdef()->countSpecs());

    // remote change caused by diffusion
    if (hostRank != myRank && !local_change) {
        if (inc <= 0) {
            std::ostringstream os;
            os << "Try to change molecule " << lidx << " by " << inc << "\n";
            os << "Fails because molecule change of receiving end should always be "
                  "non-negative.\n";
            ProgErrLog(os.str());
        }

        bufferLocations[lidx] = pRDEF->registerRemoteMoleculeChange_(
            hostRank, bufferLocations[lidx], SUB_TET, pIdx.get(), lidx, inc);
        // does not need to check sync
    }

    // local change
    else {
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

void TetRDEF::setClamped(solver::spec_local_id lidx, bool clamp) {
    if (clamp == true) {
        pPoolFlags[lidx] |= CLAMPED;
    } else {
        pPoolFlags[lidx] &= ~CLAMPED;
    }
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setOverlap(double overlap) {
    AssertLog(overlap >= 0.0);

    // Need to catch very tiny double differences
    if (overlap - pVol > double(0)) {
        AssertLog(overlap - pVol < std::numeric_limits<double>::epsilon());
        pOverlap = pVol;
    } else if (overlap / pVol < std::numeric_limits<double>::epsilon()) {
        pOverlap = 0.0;
    } else {
        pOverlap = overlap;
    }
}

////////////////////////////////////////////////////////////////////////////////

solver::linkspec_individual_id TetRDEF::getNextLinkSpecUniqueIndex() const {
    return solverRDEF()->getNextLinkspecIndividualID_();
}

////////////////////////////////////////////////////////////////////////////////

Reac& TetRDEF::reac(solver::reac_local_id lidx) const {
    AssertLog(lidx < compdef()->countReacs());
    return *dynamic_cast<Reac*>(pKProcs[lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

Diff& TetRDEF::diff(solver::diff_local_id lidx) const {
    AssertLog(lidx < compdef()->countDiffs());
    return *dynamic_cast<Diff*>(pKProcs[compdef()->countReacs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

VesBind& TetRDEF::vesbind(solver::vesbind_local_id lidx) const {
    AssertLog(lidx < compdef()->countVesBinds());
    return *dynamic_cast<VesBind*>(
        pKProcs[compdef()->countReacs() + compdef()->countDiffs() + lidx.get()]);
}

////////////////////////////////////////////////////////////////////////////////

int TetRDEF::getTetDirection(tetrahedron_global_id tidx) {
    for (uint i = 0; i < 4; i++) {
        if (pTets[i] == tidx) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::addNewLinkedSpecs(solver::linkspec_global_id linkspec1_global_id,
                                solver::linkspec_global_id linkspec2_global_id,
                                solver::linkspec_individual_id linkspec1_unique_id,
                                solver::linkspec_individual_id linkspec2_unique_id,
                                VesProxy* ves1,
                                VesProxy* ves2,
                                math::position_abs linkspec1_pos_abs,
                                math::position_abs linkspec2_pos_abs,
                                double min_length,
                                double max_length) {
    struct LinkSpec_newpair ls_pair;

    AssertLog(max_length > min_length);
    AssertLog(min_length >= 0.0);
    AssertLog(linkspec1_unique_id != linkspec2_unique_id);
    AssertLog(ves1 != ves2);

    // Do any more checks??

    ls_pair.linkspec1_global_id = linkspec1_global_id;
    ls_pair.linkspec2_global_id = linkspec2_global_id;

    ls_pair.linkspec1_individual_id = linkspec1_unique_id;
    ls_pair.linkspec2_individual_id = linkspec2_unique_id;

    ls_pair.vesicle1_individual_id = ves1->getUniqueIndex();
    ls_pair.vesicle2_individual_id = ves2->getUniqueIndex();

    ls_pair.linkspec1_pos_absolute = linkspec1_pos_abs;
    ls_pair.linkspec2_pos_absolute = linkspec2_pos_abs;

    ls_pair.min_length = min_length;
    ls_pair.max_length = max_length;

    pLinkedSpecsUpd.emplace_back(ls_pair);
}

////////////////////////////////////////////////////////////////////////////////
// MPISTEPS
bool TetRDEF::getInHost() const {
    return hostRank == myRank;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setHost(int host, int rank) {
    hostRank = host;
    myRank = rank;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setSolverRDEF(TetVesicleRDEF* solver) {
    pRDEF = solver;
}

////////////////////////////////////////////////////////////////////////////////

TetVesicleRDEF* TetRDEF::solverRDEF() const {
    return pRDEF;
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::resetPoolOccupancy() {
    std::fill(pPoolOccupancy.begin(), pPoolOccupancy.end(), 0.0);
    std::fill(pLastUpdate.begin(), pLastUpdate.end(), 0.0);
}

////////////////////////////////////////////////////////////////////////////////

void TetRDEF::setupBufferLocations() {
    uint nspecs = pCompdef->countSpecs();
    bufferLocations.container().assign(nspecs, solver::spec_local_id::unknown_value());
}

}  // namespace steps::mpi::tetvesicle
