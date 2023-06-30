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

#include "mpi/tetvesicle/sdiff.hpp"

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/mpi_common.hpp"
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/tetvesicle_rdef.hpp"
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "solver/diffdef.hpp"
#include "solver/patchdef.hpp"
#include "util/checkpointing.hpp"
#include "util/common.hpp"

// third party headers
#include "easylogging++.h"

namespace steps::mpi::tetvesicle {

SDiff::SDiff(solver::SurfDiffdef* sdef, TriRDEF* tri)
    : pSDiffdef(sdef)
    , pTri(tri)
    , pDcst(0.0)
    , pScaledDcst(0.0) {
    AssertLog(pSDiffdef != nullptr);
    AssertLog(pTri != nullptr);
    pType = KP_SDIFF;

    std::array<TriRDEF*, 3> next = {
        pTri->nextTri(0),
        pTri->nextTri(1),
        pTri->nextTri(2),
    };

    ligGIdx = pSDiffdef->lig();
    solver::Patchdef* pdef = pTri->patchdef();
    lidxTri = pdef->specG2L(ligGIdx);

    // Precalculate part of the scaled diffusion constant.
    solver::surfdiff_local_id ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);
    pDcst = dcst;

    double d[3] = {0.0, 0.0, 0.0};

    for (uint i = 0; i < 3; ++i) {
        pSDiffBndDirection[i] = pTri->getSDiffBndDirection(i);
        if (next[i] != nullptr) {
            pNeighbPatchLidx[i] = next[i]->patchdef()->specG2L(ligGIdx);
        }

        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if (!pSDiffBndDirection[i] && next[i]->patchdef() == pdef) {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
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
        pNonCDFSelector = {
            d[0] / pScaledDcst,
            d[1] / pScaledDcst,
            d[2] / pScaledDcst,
        };
    }
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, directionalDcsts);
    util::checkpoint(cp_file, pScaledDcst);
    util::checkpoint(cp_file, pDcst);
    util::checkpoint(cp_file, pNonCDFSelector);
    util::checkpoint(cp_file, pSDiffBndActive);
    util::checkpoint(cp_file, pSDiffBndDirection);
    util::checkpoint(cp_file, pNeighbPatchLidx);
    util::checkpoint(cp_file, pDirections);
    util::checkpoint(cp_file, pNdirections);
    KProc::checkpoint(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::restore(std::fstream& cp_file) {
    util::restore(cp_file, directionalDcsts);
    util::restore(cp_file, pScaledDcst);
    util::restore(cp_file, pDcst);
    util::restore(cp_file, pNonCDFSelector);
    util::restore(cp_file, pSDiffBndActive);
    util::compare(cp_file, pSDiffBndDirection);
    util::compare(cp_file, pNeighbPatchLidx);
    util::restore(cp_file, pDirections);
    util::restore(cp_file, pNdirections);
    KProc::restore(cp_file);
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::setupDeps() {
    // We will check all KProcs of the following simulation elements:
    //   * the 'source' triangle
    //   * any neighbouring tetrahedrons -- WHY??
    //
    // But also in the possible 'destination' triangles (leading to
    // four different dependency lists, each containing a copy of the
    // dependencies in the 'source' tet):
    //   * any neighbouring triangles
    //   * any neighbouring tetrahedrons of these neighbouring tris -- WHY??
    //

    AssertLog(pTri->getInHost());

    std::set<solver::kproc_global_id> remote;
    std::set<solver::kproc_global_id> remote_all;

    KProcPSet local;
    KProcPSet local_all;

    // Search for dependencies in the 'source' triangle.

    uint nkprocs = pTri->countKProcs();
    for (uint sk = 0; sk < nkprocs; sk++) {
        // Check locally.
        if (pTri->KProcDepSpecTri(sk, pTri, ligGIdx)) {
            local.insert(pTri->getKProc(sk));
        }
    }

    // Check the neighbouring tetrahedrons.
    {
        TetRDEF* itet = pTri->iTet();
        if (itet != nullptr) {
            if (pTri->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron "
                   << itet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = itet->countKProcs();
            for (uint k = 0; k < nkprocs; k++) {
                if (itet->KProcDepSpecTri(k, pTri, ligGIdx)) {
                    local.insert(itet->getKProc(k));
                }
            }
        }
    }
    {
        TetRDEF* otet = pTri->oTet();
        if (otet != nullptr) {
            if (pTri->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron "
                   << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            nkprocs = otet->countKProcs();
            for (uint k = 0; k < nkprocs; k++) {
                if (otet->KProcDepSpecTri(k, pTri, ligGIdx)) {
                    local.insert(otet->getKProc(k));
                }
            }
        }
    }

    // Search for dependencies in neighbouring triangles.
    for (uint i = 0; i < 3; ++i) {
        // Fetch next triangle, if it exists.
        TriRDEF* next = pTri->nextTri(i);
        if (next == nullptr) {
            continue;
        }

        if (next->getHost() != pTri->getHost()) {
            if (pSDiffBndDirection[i] || next->patchdef() == pTri->patchdef()) {
                pTri->solverRDEF()->addNeighHost_(next->getHost());
                pTri->solverRDEF()->registerBoundaryTri_(next);
            }
        }

        // Copy local dependencies.
        KProcPSet local2(local.begin(), local.end());
        std::set<solver::kproc_global_id> remote2;

        // Find the ones 'locally' in the next tri.
        nkprocs = next->countKProcs();
        auto startKProcIdx = next->getStartKProcIdx();
        for (uint sk = 0; sk < nkprocs; sk++) {
            // Check locally.
            if (next->KProcDepSpecTri(sk, next, ligGIdx)) {
                if (next->getHost() == pTri->getHost()) {
                    local2.insert(next->getKProc(sk));
                } else {
                    remote2.emplace(startKProcIdx + sk);
                }
            }
        }

        // Fetch inner tetrahedron, if it exists.
        TetRDEF* itet = next->iTet();
        if (itet != nullptr) {
            if (next->getHost() != itet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << pTri->idx() << " and its compartment tetrahedron "
                   << itet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }

            // Find deps.
            nkprocs = itet->countKProcs();
            startKProcIdx = itet->getStartKProcIdx();
            for (uint k = 0; k < nkprocs; k++) {
                if (itet->KProcDepSpecTri(k, next, ligGIdx)) {
                    if (itet->getHost() == pTri->getHost()) {
                        local2.insert(itet->getKProc(k));
                    } else {
                        remote2.emplace(startKProcIdx + k);
                    }
                }
            }
        }

        // Fetch outer tetrahedron, if it exists.
        TetRDEF* otet = next->oTet();
        if (otet != nullptr) {
            if (next->getHost() != otet->getHost()) {
                std::ostringstream os;
                os << "Patch triangle " << next->idx() << " and its compartment tetrahedron "
                   << otet->idx() << " belong to different hosts.\n";
                NotImplErrLog(os.str());
            }
            // Find deps.
            nkprocs = otet->countKProcs();
            startKProcIdx = otet->getStartKProcIdx();
            for (uint k = 0; k < nkprocs; k++) {
                if (otet->KProcDepSpecTri(k, next, ligGIdx)) {
                    if (otet->getHost() == pTri->getHost()) {
                        local2.insert(otet->getKProc(k));
                    } else {
                        remote2.emplace(startKProcIdx + k);
                    }
                }
            }
        }

        localUpdVec[i].assign(local2.begin(), local2.end());
        remoteUpdVec[i].assign(remote2.begin(), remote2.end());

        local_all.insert(local2.begin(), local2.end());
        remote_all.insert(remote2.begin(), remote2.end());
    }

    localAllUpdVec.assign(local_all.begin(), local_all.end());
    remoteAllUpdVec.assign(remote_all.begin(), remote_all.end());

    // pSpecChange.insert(ligGIdx);
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::reset() {
    resetExtent();

    pSDiffBndActive = {false, false, false};

    solver::surfdiff_local_id ldidx = pTri->patchdef()->surfdiffG2L(pSDiffdef->gidx());
    double dcst = pTri->patchdef()->dcst(ldidx);

    setDcst(dcst);

    setActive(true);

    crData.recorded = false;
    crData.pow = 0;
    // cannot reset because the position is used in solver
    // crData.pos = 0;
    crData.rate = 0.0;
}

////////////////////////////////////////////////////////////////////////////////

double SDiff::dcst(int direction) {
    if (directionalDcsts.find(direction) != directionalDcsts.end()) {
        return directionalDcsts[direction];
    } else {
        return pDcst;
    }
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::setDcst(double dcst) {
    AssertLog(dcst >= 0.0);
    pDcst = dcst;
    directionalDcsts.clear();

    std::array<TriRDEF*, 3> next{pTri->nextTri(0), pTri->nextTri(1), pTri->nextTri(2)};

    // Reset this stuff- may have been created before, may not have been (if
    // original dcst was 0)
    pNdirections = 0;
    pDirections.clear();

    std::array<double, 3> d{0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    for (uint i = 0; i < d.size(); ++i) {
        // Compute the scaled diffusion constant.
        // Need to here check if the direction is a diffusion boundary
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) ||
                (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef())) {
                d[i] = (pTri->length(i) * dcst) / (pTri->area() * dist);
            }
        }
        if (d[i] > 0.0) {
            pScaledDcst += d[i];
            pDirections.push_back(i);
            pNdirections += 1;
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pNonCDFSelector = {0.0, 0.0, 0.0};
    } else {
        pNonCDFSelector = {d[0] / pScaledDcst, d[1] / pScaledDcst, d[2] / pScaledDcst};
    }
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::setDirectionDcst(int direction, double dcst) {
    AssertLog(direction < 3);
    AssertLog(direction >= 0);
    AssertLog(dcst >= 0.0);
    directionalDcsts[direction] = dcst;

    TriRDEF* next[3] = {pTri->nextTri(0), pTri->nextTri(1), pTri->nextTri(2)};

    // Automatically activate boundary diffusion if necessary TODO check this
    if (pSDiffBndDirection[direction] == true) {
        pSDiffBndActive[direction] = true;
    }

    double d[3] = {0.0, 0.0, 0.0};
    pScaledDcst = 0.0;

    pNdirections = 0;
    pDirections.clear();

    for (uint i = 0; i < 3; ++i) {
        // if directional diffusion dcst exists use directional dcst, else use the
        // standard one
        double dist = pTri->dist(i);
        if ((dist > 0.0) && (next[i] != nullptr)) {
            if ((pSDiffBndDirection[i] && pSDiffBndActive[i]) ||
                (!pSDiffBndDirection[i] && next[i]->patchdef() == pTri->patchdef())) {
                auto search_result = directionalDcsts.find(i);
                if (search_result != directionalDcsts.end()) {
                    d[i] = (pTri->length(i) * search_result->second) / (pTri->area() * dist);
                } else {
                    // This part must use the default pDcst, not the function argument
                    // dcst
                    d[i] = (pTri->length(i) * pDcst) / (pTri->area() * dist);
                }
            }
        }
        if (d[i] > 0.0) {
            pScaledDcst += d[i];
            pDirections.push_back(i);
            pNdirections += 1;
        }
    }

    // Should not be negative!
    AssertLog(pScaledDcst >= 0);

    // Setup the selector distribution.
    if (pScaledDcst == 0.0) {
        pNonCDFSelector = {0.0, 0.0, 0.0};
    } else {
        pNonCDFSelector = {d[0] / pScaledDcst, d[1] / pScaledDcst, d[2] / pScaledDcst};
    }
}

////////////////////////////////////////////////////////////////////////////////

double SDiff::rate(TetVesicleRDEF* /*solver*/) {
    if (inactive()) {
        return 0.0;
    }

    // Compute the rate.
    double rate = pScaledDcst * static_cast<double>(pTri->pools()[lidxTri]);
    AssertLog(std::isnan(rate) == false);

    // Return.
    return rate;
}

////////////////////////////////////////////////////////////////////////////////

int SDiff::apply(const rng::RNGptr& rng) {
    // Apply local change.
    uint local = pTri->pools()[lidxTri];
    bool clamped = pTri->clamped(lidxTri);

    // pConnectedTets[0] = pTri->iTet, pConnectedTets[1] = pTri->oTet
    // rest (re)set to nullptr
    // std::fill (pConnectedTets.begin()+2, pConnectedTets.end(), nullptr);

    if (clamped == false) {
        if (local == 0) {
            return -2;
        }  // no molecule left, no diffusion
    }

    // We should have a direction ergo pConnected tets size
    // should be bigger than just 2, the source tri connected tets. Use
    // pConnectedTets[2, 3] below so do this assert here
    // AssertLog(pConnectedTets.size() >= 4);

    // Apply change in next voxel: select a direction.
    double sel = rng->getUnfEE();

    int iSel = 0;
    double CDFSelector = 0.0;
    for (; iSel < 2; ++iSel) {
        CDFSelector += pNonCDFSelector[iSel];
        if (sel < CDFSelector) {
            break;
        }
    }

    // Direction iSel.
    TriRDEF* nexttri = pTri->nextTri(iSel);

    AssertLog(nexttri != nullptr);
    AssertLog(pNeighbPatchLidx[iSel].valid());

    if (nexttri->clamped(pNeighbPatchLidx[iSel]) == false) {
        nexttri->incCount(pNeighbPatchLidx[iSel], 1);
    }

    if (clamped == false) {
        pTri->incCount(lidxTri, -1);
    }

    rExtent++;

    return iSel;
}

///////////////////////////////////////////////////////////////////////////////

int SDiff::apply(const rng::RNGptr& rng, uint nmolcs) {
    // Apply local change.
    uint local = pTri->pools()[lidxTri];
    bool clamped = pTri->clamped(lidxTri);

    if (clamped == false) {
        if (local == 0) {
            return -2;
        }
    }

    AssertLog(pNdirections >= 1);

    // Multinomial
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
            TriRDEF* nexttri = pTri->nextTri(direction);

            AssertLog(nexttri != nullptr);
            AssertLog(pNeighbPatchLidx[direction].valid());

            if (nexttri->clamped(pNeighbPatchLidx[direction]) == false) {
                // 2 connected tets for source, 2 per direction,
                // hence these indices
                // pConnectedTets[2+i*2] = nexttri->oTet();
                // pConnectedTets[3+i*2] = nexttri->iTet();
                nexttri->incCount(pNeighbPatchLidx[direction], molcsthisdir);
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
        TriRDEF* nexttri = pTri->nextTri(direction);

        AssertLog(nexttri != nullptr);
        AssertLog(pNeighbPatchLidx[direction].valid());

        if (nexttri->clamped(pNeighbPatchLidx[direction]) == false) {
            // pConnectedTets[pNdirections*2] = nexttri->oTet();
            // pConnectedTets[pNdirections*2+1] = nexttri->iTet();
            nexttri->incCount(pNeighbPatchLidx[direction], molcsthisdir);
        }

        molcs_moved += molcsthisdir;
    }

    AssertLog(molcs_moved == nmolcs);

    if (clamped == false) {
        pTri->incCount(lidxTri, -nmolcs);
    }

    rExtent += nmolcs;

    return -1;
}

////////////////////////////////////////////////////////////////////////////////

std::vector<solver::kproc_global_id> const& SDiff::getRemoteUpdVec(int direction) const {
    if (direction == -1) {
        return remoteAllUpdVec;
    } else if (direction == -2) {
        return idxEmptyvec;
    } else {
        return remoteUpdVec[direction];
    }
}

////////////////////////////////////////////////////////////////////////////////

std::vector<KProc*> const& SDiff::getLocalUpdVec(int direction) const {
    if (direction == -1) {
        return localAllUpdVec;
    } else if (direction == -2) {
        return pEmptyvec;
    } else {
        return localUpdVec[direction];
    }
}

////////////////////////////////////////////////////////////////////////////////

void SDiff::setSDiffBndActive(uint i, bool active) {
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    // Only need to update if the flags are changing
    if (pSDiffBndActive[i] != active) {
        pSDiffBndActive[i] = active;
        setDcst(pDcst);
    }
}

////////////////////////////////////////////////////////////////////////////////

bool SDiff::getSDiffBndActive(uint i) const {
    AssertLog(i < 3);
    AssertLog(pSDiffBndDirection[i] == true);

    return pSDiffBndActive[i];
}

}  // namespace steps::mpi::tetvesicle
