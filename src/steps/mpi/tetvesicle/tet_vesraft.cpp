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

#include "mpi/tetvesicle/tet_vesraft.hpp"

// Standard library & STL headers.
#include <algorithm>
#include <cmath>
#include <functional>

// STEPS headers.
#include "mpi/tetvesicle/exocytosis.hpp"

#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "solver/compdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

TetVesRaft::TetVesRaft(tetrahedron_global_id idx,
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
                       math::point3d baryc)
    : pIdx(idx)
    , pCompdef(cdef)
    , pVol(vol)
    , pOverlap(0.0)
    , pReducedVol(vol)
    , pPosition(baryc)
    , pNextTet()
    , pAreas()
    , pDist() {
    AssertLog(pCompdef != nullptr);
    AssertLog(pVol > 0.0);

    AssertLog(a0 > 0.0 && a1 > 0.0 && a2 > 0.0 && a3 > 0.0);
    AssertLog(d0 >= 0.0 && d1 >= 0.0 && d2 >= 0.0 && d3 >= 0.0);

    // Based on compartment definition, build other structures.
    uint nspecs = compdef()->countSpecs();
    pPoolCount.container().resize(nspecs);

    // At this point we don't have neighbouring tet pointers,
    // but we can store their indices
    for (uint i = 0; i <= pNextTet.size(); ++i) {
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
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolCount);
    // NOTE no checkpoint for pOverlap, pReducedVol because this comes when vesicles are restored
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::setNextTet(uint i, TetVesRaft* t) {
    // Now adding all tets, even those from other compartments, due to the
    // diffusion boundaries
    pNextTet[i] = t;
    pNextTris[i] = nullptr;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::setNextTri(uint i, TriVesRaft* t) {
    AssertLog(pNextTris.size() == 4);
    AssertLog(i <= 3);

    pNextTet[i] = nullptr;
    pNextTris[i] = t;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::reset() {
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);

    pOverlap = 0.0;
    pReducedVol = pVol;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::setCount(solver::spec_local_id lidx, uint count) {
    // Count has changed, need to correct pool factor

    // This function is used for updates that do not require remote sync,
    // such as user setCount or reaction

    pPoolCount.at(lidx) = count;
}

////////////////////////////////////////////////////////////////////////////////

void TetVesRaft::changeOverlap(double overlap) {
    double new_overlap = pOverlap + overlap;

    // Need to catch very tiny double differences, full overlap
    if (new_overlap - pVol > double(0)) {
        AssertLog(new_overlap - pVol < std::numeric_limits<double>::epsilon());
        pOverlap = pVol;
        pReducedVol = solver::TINY_VOLUME;
    }
    // Need to catch very tiny double differences, no overlap
    else if (new_overlap / pVol < std::numeric_limits<double>::epsilon()) {
        pOverlap = 0.0;
        pReducedVol = pVol;
    } else {
        pOverlap = new_overlap;
        pReducedVol = pVol - pOverlap;
    }
}

////////////////////////////////////////////////////////////////////////////////

int TetVesRaft::getTetDirection(tetrahedron_global_id tidx) {
    for (uint i = 0; i < pTets.size(); i++) {
        if (pTets[i] == tidx) {
            return static_cast<int>(i);
        }
    }
    return -1;
}

}  // namespace steps::mpi::tetvesicle
