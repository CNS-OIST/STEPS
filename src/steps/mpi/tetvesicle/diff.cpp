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

#include "mpi/tetvesicle/diff.hpp"

// Standard library & STL headers.
#include <iostream>
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tet_rdef.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "solver/compdef.hpp"
#include "solver/diffdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

Diff::Diff(solver::Diffdef* ddef, TetRDEF* tet)
    : pDiffdef(ddef)
    , pTet(tet)
    , pDcst(0.0)
    , pScaledDcst(0.0) {
    AssertLog(pDiffdef != nullptr);
    AssertLog(pTet != nullptr);
    pType = KP_DIFF;

    const auto& next = pTet->nextTets();

    ligGIdx = pDiffdef->lig();
    solver::Compdef* cdef = pTet->compdef();
    lidxTet = cdef->specG2L(ligGIdx);

    // Precalculate part of the scaled diffusion constant.
    solver::diff_local_id ldidx(pTet->compdef()->diffG2L(pDiffdef->gidx()));
    double dcst = pTet->compdef()->dcst(ldidx);
    pDcst = dcst;

    double d[4] = {0.0, 0.0, 0.0, 0.0};

    for (uint i = 0; i < 4; ++i) {
        pDiffBndDirection[i] = pTet->getDiffBndDirection(i);
        if (next[i] != nullptr) {
            pNeighbCompLidx[i] = next[i]->compdef()->specG2L(ligGIdx);
        }

        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if (!pDiffBndDirection[i] && next[i]->compdef() == cdef) {
                d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
                if (d[i] > 0.0) {
                    pScaledDcst += d[i];
                    pDirections.push_back(i);
                    pNdirections += 1;
                }
            }
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst > 0.0) {
        pNonCDFSelector = {d[0] / pScaledDcst,
                           d[1] / pScaledDcst,
                           d[2] / pScaledDcst,
                           d[3] / pScaledDcst};
    }
}

////////////////////////////////////////////////////////////////////////////////

Diff::~Diff() = default;

////////////////////////////////////////////////////////////////////////////////

void Diff::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, directionalDcsts);
    util::checkpoint(cp_file, pScaledDcst);
    util::checkpoint(cp_file, pDcst);
    util::checkpoint(cp_file, pNonCDFSelector);
    util::checkpoint(cp_file, pDiffBndActive);
    util::checkpoint(cp_file, pDiffBndDirection);
    util::checkpoint(cp_file, pNeighbCompLidx);
    util::checkpoint(cp_file, pDirections);
    util::checkpoint(cp_file, pNdirections);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Diff::restore(std::fstream& cp_file) {
    util::restore(cp_file, directionalDcsts);
    util::restore(cp_file, pScaledDcst);
    util::restore(cp_file, pDcst);
    util::restore(cp_file, pNonCDFSelector);
    util::restore(cp_file, pDiffBndActive);
    util::compare(cp_file, pDiffBndDirection);
    util::compare(cp_file, pNeighbCompLidx);
    util::restore(cp_file, pDirections);
    util::restore(cp_file, pNdirections);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setupDeps() {
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' tetrahedron
    //   * any neighbouring triangles
    //
    // But also in the possible 'destination' tetrahedrons (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring tetrahedrons
    //   * any neighbouring triangles of these neighbouring tets
    //
    // Since there can be no diffusion between tetrahedrons blocked by
    // a triangle, there is no need to filter out duplicate dependent
    // kprocs.

    // Search for dependencies in the 'source' tetrahedron.

    AssertLog(pTet->getInHost());

    std::set<solver::kproc_global_id> remote;
    std::set<solver::kproc_global_id> remote_all;

    KProcPSet local;
    KProcPSet local_all;

    uint nkprocs = pTet->countKProcs();
    for (uint k = 0; k < nkprocs; k++) {
        if (pTet->KProcDepSpecTet(k, pTet, ligGIdx)) {
            local.insert(pTet->getKProc(k));
        }
    }
    // Check the neighbouring triangles.
    for (uint i = 0; i < 4; ++i) {
        TriRDEF* next = pTet->nextTri(i);
        if (next == nullptr) {
            continue;
        }

        // next tri has to be in the same host to prevent
        // cross process surface reaction
        if (next->getHost() != pTet->getHost()) {
            std::ostringstream os;
            os << "Patch triangle " << next->idx() << " and its compartment ";
            os << "tetrahedron " << pTet->idx() << " belong to different hosts.\n";
            NotImplErrLog(os.str());
        }

        nkprocs = next->countKProcs();
        for (uint sk = 0; sk < nkprocs; sk++) {
            if (next->KProcDepSpecTet(sk, pTet, ligGIdx)) {
                local.insert(next->getKProc(sk));
            }
        }
    }

    // Search for dependencies in neighbouring tetrahedra.
    for (uint i = 0; i < 4; ++i) {
        // Fetch next tetrahedron, if it exists.
        TetRDEF* next = pTet->nextTet(i);
        if (next == nullptr) {
            continue;
        }
        if (pTet->nextTri(i) != nullptr) {
            continue;
        }

        // Copy local dependencies.
        KProcPSet local2(local.begin(), local.end());
        std::set<solver::kproc_global_id> remote2;

        // Find the ones 'locally' in the next tet.
        nkprocs = next->countKProcs();
        auto startKProcIdx = next->getStartKProcIdx();

        if (next->getHost() != pTet->getHost()) {
            if (pDiffBndDirection[i] || next->compdef() == pTet->compdef()) {
                pTet->solverRDEF()->addNeighHost_(next->getHost());
                pTet->solverRDEF()->registerBoundaryTet_(next);
            }
        }

        for (uint k = 0; k < nkprocs; k++) {
            if (next->KProcDepSpecTet(k, next, ligGIdx)) {
                // if next tet has same host as pTet, store dependent kp as pointer
                if (next->getHost() == pTet->getHost()) {
                    local2.insert(next->getKProc(k));
                }
                // if not store as index
                else {
                    remote2.emplace(startKProcIdx + k);
                }
            }
        }

        // Find deps in neighbouring triangles in the next tet.
        // As said before, this cannot logically include the shared
        // triangle.
        for (uint j = 0; j < 4; ++j) {
            // Fetch next triangle, if it exists.
            TriRDEF* next2 = next->nextTri(j);
            if (next2 == nullptr) {
                continue;
            }

            if (next2->getHost() != next->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next2->idx() << " and its compartment";
                os << " tetrahedron " << next->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            // Find deps.
            nkprocs = next2->countKProcs();
            startKProcIdx = next2->getStartKProcIdx();
            for (uint sk = 0; sk < nkprocs; sk++) {
                if (next2->KProcDepSpecTet(sk, next, ligGIdx)) {
                    if (next2->getHost() == pTet->getHost()) {
                        local2.insert(next2->getKProc(sk));
                    } else {
                        remote2.emplace(startKProcIdx + sk);
                    }
                }
            }
        }

        // Copy the set to the update vector.
        localUpdVec[i].assign(local2.begin(), local2.end());
        remoteUpdVec[i].assign(remote2.begin(), remote2.end());
        local_all.insert(local2.begin(), local2.end());
        remote_all.insert(remote2.begin(), remote2.end());
    }

    localAllUpdVec.assign(local_all.begin(), local_all.end());
    remoteAllUpdVec.assign(remote_all.begin(), remote_all.end());
}

////////////////////////////////////////////////////////////////////////////////

void Diff::reset() {
    resetExtent();

    // NOTE: These must become the dcst calculation for obvious reasons
    pDiffBndActive = {false, false, false, false};

    solver::diff_local_id ldidx = pTet->compdef()->diffG2L(pDiffdef->gidx());
    double dcst = pTet->compdef()->dcst(ldidx);

    // directional dcst will also be clear by setDcst
    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    // pos cannot be reset because their positions in pDiffs will be used to swap
    // crData.pos = 0;
    crData.rate = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDiffBndActive(uint i, bool active) {
    AssertLog(i < 4);
    AssertLog(pDiffBndDirection[i] == true);

    // Only need to update if the flags are changing
    if (pDiffBndActive[i] != active) {
        pDiffBndActive[i] = active;
        setDcst(pDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool Diff::getDiffBndActive(uint i) const {
    AssertLog(i < 4);
    AssertLog(pDiffBndDirection[i] == true);

    return pDiffBndActive[i];
}

////////////////////////////////////////////////////////////////////////////////

double Diff::dcst(int direction) {
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    } else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDcst(double dcst) {
    AssertLog(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    TetRDEF* next[4] = {pTet->nextTet(0), pTet->nextTet(1), pTet->nextTet(2), pTet->nextTet(3)};

    // Reset this stuff- may have been created before, may not have been (if
    // original dcst was 0)
    pNdirections = 0;
    pDirections.clear();

    double d[4] = {0.0, 0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    for (uint i = 0; i < 4; ++i) {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            // Don't diffuse to tets of zero volume
            if (next[i]->vol() == solver::TINY_VOLUME) {
                continue;
            }
            if ((pDiffBndDirection[i] && pDiffBndActive[i]) ||
                (!pDiffBndDirection[i] && next[i]->compdef() == pTet->compdef())) {
                d[i] = (pTet->area(i) * dcst) / (pTet->vol() * dist);
            }
        }
        if (d[i] > 0.0) {
            pScaledDcst += d[i];
            pDirections.push_back(i);
            pNdirections += 1;
        }
    }

    // pConnectedTets.resize(1+pNdirections);
    // pConnectedTets[0] = pTet;

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pNonCDFSelector = {0.0, 0.0, 0.0, 0.0};
    } else {
        pNonCDFSelector = {d[0] / pScaledDcst,
                           d[1] / pScaledDcst,
                           d[2] / pScaledDcst,
                           d[3] / pScaledDcst};
    }
}

////////////////////////////////////////////////////////////////////////////////

void Diff::setDirectionDcst(int direction, double dcst) {
    AssertLog(direction < 4);
    AssertLog(direction >= 0);
    AssertLog(dcst >= 0.0);
    directionalDcsts[direction] = dcst;

    recalcDcst();
}

////////////////////////////////////////////////////////////////////////////////

void Diff::recalcDcst() {
    TetRDEF* next[4] = {pTet->nextTet(0), pTet->nextTet(1), pTet->nextTet(2), pTet->nextTet(3)};

    // Reset this stuff- may have been created before, may not have been (if
    // original dcst was 0)
    pNdirections = 0;
    pDirections.clear();

    double d[4] = {0.0, 0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    for (uint i = 0; i < 4; ++i) {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTet->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            // Don't diffuse to tets of zero volume
            if (next[i]->vol() == solver::TINY_VOLUME) {
                continue;
            }

            if ((pDiffBndDirection[i] && pDiffBndActive[i]) ||
                (!pDiffBndDirection[i] && next[i]->compdef() == pTet->compdef())) {
                if (directionalDcsts.find(i) != directionalDcsts.end()) {
                    d[i] = (pTet->area(i) * directionalDcsts[i]) / (pTet->vol() * dist);
                } else {
                    d[i] = (pTet->area(i) * pDcst) / (pTet->vol() * dist);
                }
            }
        }
        if (d[i] > 0.0) {
            pScaledDcst += d[i];
            pDirections.push_back(i);
            pNdirections += 1;
        }
    }

    // pConnectedTets.resize(1+pNdirections);
    // pConnectedTets[0] = pTet;

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pNonCDFSelector = {0.0, 0.0, 0.0, 0.0};

    } else {
        pNonCDFSelector = {d[0] / pScaledDcst,
                           d[1] / pScaledDcst,
                           d[2] / pScaledDcst,
                           d[3] / pScaledDcst};
    }
}

////////////////////////////////////////////////////////////////////////////////

double Diff::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // This is necessary for the hybrid solver- first the rates must be
    // recalculated because volumes may have changed by vesicle overlap.
    recalcDcst();

    // Little time-saver. Huge for sparsely populated meshes.
    uint pool = pTet->pools()[lidxTet];
    if (pool == 0) {
        return 0.0;
    }

    // Compute the rate.
    double rate = (pScaledDcst) * static_cast<double>(pool);
    AssertLog(std::isnan(rate) == false);

    return rate;
}

////////////////////////////////////////////////////////////////////////////////

int Diff::apply(const rng::RNGptr& rng) {
    // Apply local change.
    bool clamped = pTet->clamped(lidxTet);

    // Remember pConnectedTets[0] = pTet, rest reset to nullptr
    // std::fill (pConnectedTets.begin()+1, pConnectedTets.end(), nullptr);

    if (clamped == false) {
        auto local = pTet->pools()[lidxTet];
        if (local == 0) {
            return -2;
        }  // no molecule left, no diffusion
    }

    // We should have a direction ergo pConnected tets size
    // should be bigger than just 1, the source tet. Use
    // pConnectedTets[1] below so do this assert here
    // AssertLog(pConnectedTets.size() > 1);

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    int iSel = 0;
    double CDFSelector = 0.0;
    for (; iSel < 3; ++iSel) {
        CDFSelector += pNonCDFSelector[iSel];
        if (sel < CDFSelector) {
            break;
        }
    }

    // Direction iSel.
    TetRDEF* nexttet = pTet->nextTet(iSel);

    AssertLog(nexttet != nullptr);
    AssertLog(pNeighbCompLidx[iSel].valid());

    if (nexttet->clamped(pNeighbCompLidx[iSel]) == false) {
        // For the vesicle kprocs, useful to know which tets were affected
        // pConnectedTets[1] = nexttet; // TODO check set this only if moved
        nexttet->incCount(pNeighbCompLidx[iSel], 1);
    }
    // else { pConnectedTets[1] = nullptr; }

    if (clamped == false) {
        pTet->incCount(lidxTet, -1);
    }

    rExtent++;

    return iSel;
}

///////////////////////////////////////////////////////////////////////////////

int Diff::apply(const rng::RNGptr& rng, uint nmolcs) {
    // Apply local change.
    bool clamped = pTet->clamped(lidxTet);

    if (clamped == false) {
        auto local = pTet->pools()[lidxTet];
        if (local == 0) {
            return -2;
        }
    }

    AssertLog(pNdirections >= 1);
    // AssertLog(pConnectedTets.size() == pNdirections+1);

    // Start with nullptr because previous data may be held and we may
    // may not visit every direction. Remember pConnectedTets[0] = pTet
    // std::fill (pConnectedTets.begin()+1, pConnectedTets.end(), nullptr);

    // Multinomial by stl
    uint molcs_moved = 0;
    for (uint i = 0; i < pNdirections - 1; ++i) {
        uint direction = pDirections[i];
        double chance = 0.0;

        double sump = 0.0;
        for (uint j = 0; j < i; ++j) {
            sump += pNonCDFSelector[pDirections[j]];
        }
        chance = pNonCDFSelector[direction] / (1.0 - sump);

        unsigned int max_molcs = nmolcs - molcs_moved;

        // To avoid cycling to the largest unsigned int
        if (chance >= 1.0) {
            chance = 1.0;
        }

        uint molcsthisdir = rng->getBinom(max_molcs, chance);

        if (molcsthisdir != 0) {
            TetRDEF* nexttet = pTet->nextTet(direction);

            AssertLog(nexttet != nullptr);
            AssertLog(pNeighbCompLidx[direction].valid());

            if (nexttet->clamped(pNeighbCompLidx[direction]) == false) {
                // 0 index is the source tet, hence the i+1
                // pConnectedTets[i+1] = nexttet;
                nexttet->incCount(pNeighbCompLidx[direction], molcsthisdir);
            }

            molcs_moved += molcsthisdir;
        }
        if (molcs_moved == nmolcs) {
            break;
        }
    }
    // last direction, chance =1
    uint direction = pDirections[pNdirections - 1];
    int molcsthisdir = nmolcs - molcs_moved;
    if (molcsthisdir != 0) {
        TetRDEF* nexttet = pTet->nextTet(direction);

        AssertLog(nexttet != nullptr);
        AssertLog(pNeighbCompLidx[direction].valid());

        if (nexttet->clamped(pNeighbCompLidx[direction]) == false) {
            // pConnectedTets[pNdirections] = nexttet;
            nexttet->incCount(pNeighbCompLidx[direction], molcsthisdir);
        }

        molcs_moved += molcsthisdir;
    }

    AssertLog(molcs_moved == nmolcs);

    if (clamped == false) {
        pTet->incCount(lidxTet, -nmolcs);
    }

    rExtent += nmolcs;

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& Diff::getRemoteUpdVec(int direction) const {
    if (direction == -1) {
        return remoteAllUpdVec;
    } else if (direction == -2) {
        return idxEmptyvec;
    } else {
        return remoteUpdVec[direction];
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& Diff::getLocalUpdVec(int direction) const {
    if (direction == -1) {
        return localAllUpdVec;
    } else if (direction == -2) {
        return pEmptyvec;
    } else {
        return localUpdVec[direction];
    }
}

}  // namespace steps::mpi::tetvesicle
