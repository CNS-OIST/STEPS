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

#pragma once

// STL headers.
#include <fstream>
#include <list>
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/tri_rdef.hpp"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
// class PatchRDEF;
class Endocytosis;

////////////////////////////////////////////////////////////////////////////////

class PatchRDEF {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    PatchRDEF(solver::Patchdef* patchdef, tetmesh::Tetmesh* mesh, TetVesicleRDEF* rdef);
    ~PatchRDEF();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // Additions and setup
    ////////////////////////////////////////////////////////////////////////

    void reset() const;

    // void setupRafts();

    ////////////////////////////////////////////////////////////////////////
    // DATA, MEMBER ACCESS
    ////////////////////////////////////////////////////////////////////////

    inline solver::Patchdef* def() const noexcept {
        return pPatchdef;
    }

    inline tetmesh::Tetmesh* mesh() const noexcept {
        return pMesh;
    }

    inline double area() const noexcept {
        return pArea;
    }

    /// Return the random number generator object.
    inline const rng::RNGptr& rng() const noexcept {
        return pRNG;
    }

    ////////////////////////////////////////////////////////////////////////
    // TRIANGLES
    ////////////////////////////////////////////////////////////////////////

    /// Checks whether Tri::patchdef() corresponds to this object's
    /// PatchDef. There is no check whether the Tri object has already
    /// been added to this Patch object before (i.e. no duplicate
    /// checking).
    ///
    void addTri(TriRDEF* tri);

    inline uint countTris() const {
        return pTris.size();
    }

    inline const TriRDEFPVec& tris() const noexcept {
        return pTris;
    }

    TriRDEF* tri(triangle_local_id tri_lidx) const noexcept {
        return pTris[tri_lidx.get()];
    }

    inline triangle_global_id triidx_L_to_G(triangle_local_id lidx) noexcept {
        return pTriidcs_L_to_G[lidx];
    }
    inline triangle_local_id triidx_G_to_L(triangle_global_id gidx) noexcept {
        return pTriidcs_G_to_L[gidx];
    }

    inline TetVesicleRDEF* solverRDEF() {
        return pRDEF;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Patchdef* pPatchdef;
    double pArea;

    tetmesh::Tetmesh* pMesh;

    // Just easier to store the solver for all the raft KProc stuff
    TetVesicleRDEF* pRDEF{nullptr};

    // Store a pointer to the RNG for convenience
    rng::RNGptr pRNG;

    TriRDEFPVec pTris;

    std::map<triangle_global_id, triangle_local_id> pTriidcs_G_to_L;
    std::map<triangle_local_id, triangle_global_id> pTriidcs_L_to_G;
};

}  // namespace steps::mpi::tetvesicle
