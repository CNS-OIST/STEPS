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

#include "mpi/tetvesicle/tri_vesraft.hpp"

// Standard library & STL headers.
#include <algorithm>
#include <cmath>
#include <functional>

// STEPS headers.
#include "math/constants.hpp"
#include "mpi/tetvesicle/tet_vesraft.hpp"
#include "mpi/tetvesicle/tetvesicle_vesraft.hpp"
#include "solver/patchdef.hpp"
#include "util/checkpointing.hpp"

namespace steps::mpi::tetvesicle {

TriVesRaft::TriVesRaft(triangle_global_id idx,
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
                       const math::point3d& trinorm)
    : pIdx(idx)
    , pPatchdef(patchdef)
    , pArea(area)
    , pInnerTet(nullptr)
    , pOuterTet(nullptr)
    , pNextTri()
    , pLengths()
    , pDist()
    , pPosition(position)
    , pNorm(trinorm)
    , pPatchVesRaft(nullptr) {
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

    const uint nspecs = pPatchdef->countSpecs();
    pPoolCount.container().resize(nspecs);
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::checkpoint(std::fstream& cp_file) {
    util::checkpoint(cp_file, pPoolCount);
    // Note pAppliedRaftgens cleared at time of call to checkpoint
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::restore(std::fstream& cp_file) {
    util::restore(cp_file, pPoolCount);
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::setInnerTet(TetVesRaft* t) {
    pInnerTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::setOuterTet(TetVesRaft* t) {
    pOuterTet = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::setNextTri(uint i, TriVesRaft* t) {
    pNextTri.at(i) = t;
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::reset() {
    std::fill(pPoolCount.begin(), pPoolCount.end(), 0);
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::setCount(solver::spec_local_id lidx, uint count) {
    AssertLog(lidx < patchdef()->countSpecs());

    pPoolCount[lidx] = count;
}

////////////////////////////////////////////////////////////////////////////////

int TriVesRaft::getTriDirection(triangle_global_id tidx) {
    for (uint i = 0; i < pTris.size(); i++) {
        if (pTris[i] == tidx) {
            return i;
        }
    }
    return -1;
}

////////////////////////////////////////////////////////////////////////////////

void TriVesRaft::setSolverVesRaft(TetVesicleVesRaft* solver) {
    pVesRaft = solver;
}
////////////////////////////////////////////////////////////////////////////////

TetVesicleVesRaft* TriVesRaft::solverVesRaft() const {
    return pVesRaft;
}

}  // namespace steps::mpi::tetvesicle
