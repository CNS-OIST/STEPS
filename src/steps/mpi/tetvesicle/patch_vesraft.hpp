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
#include "mpi/tetvesicle/tri_vesraft.hpp"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations.
class TetVesicleVesRaft;
class Endocytosis;

////////////////////////////////////////////////////////////////////////////////

class PatchVesRaft {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    PatchVesRaft(solver::Patchdef* patchdef, tetmesh::Tetmesh* mesh, TetVesicleVesRaft* vesraft);
    ~PatchVesRaft();

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

    void reset();

    void setupRafts();

    void setupEndocyticZones();

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
    void addTri(TriVesRaft* tri);

    inline uint countTris() const {
        return pTris.size();
    }

    TriVesRaft* pickTriByArea(double rand01) const;

    inline const auto& tris() const noexcept {
        return pTris;
    }

    TriVesRaft* tri(triangle_local_id tri_lidx) const noexcept {
        return pTris[tri_lidx];
    }

    inline triangle_global_id triidx_L_to_G(triangle_local_id lidx) noexcept {
        return pTriidcs_L_to_G[lidx];
    }
    inline triangle_local_id triidx_G_to_L(triangle_global_id gidx) noexcept {
        return pTriidcs_G_to_L[gidx];
    }

    ////////////////////////////////////////////////////////////////////////
    // ENDOCYTOSIS and ZONES
    ////////////////////////////////////////////////////////////////////////

    Endocytosis* getEndocytosis(solver::endocytosis_local_id endoidx, const std::string& zone);

    ////////////////////////////////////////////////////////////////////////
    // RAFTS
    ////////////////////////////////////////////////////////////////////////

    std::map<solver::raft_global_id, std::list<Raft*>> const& getAllRafts() const noexcept {
        return pRafts;
    }
    // Raft can either be added at a given position (intended for exocytosis)
    // or at a random position within the patch
    // returns unknown id if Raft not added
    solver::raft_individual_id addRaft(solver::Raftdef* raftdef, TriVesRaft* tri);
    solver::raft_individual_id addRaft(solver::raft_global_id ridx, TriVesRaft* tri);

    void removeOneRaft(Raft* raft);

    // Set/get count anywhere and everywhere
    void setRaftCount(solver::raft_global_id ridx, uint count);
    uint getRaftCount(solver::raft_global_id ridx) const;

    // In/from specific triangles
    void setRaftCount(solver::raft_global_id ridx, TriVesRaft* tri, uint count);
    uint getRaftCount(solver::raft_global_id ridx, TriVesRaft* tri) const;

    std::vector<solver::raft_individual_id> getRaftIndices(solver::raft_global_id ridx) const;

    // Returns as a map of individual raft (by unique index) to the number of this
    // species present in that Raft
    std::map<solver::raft_individual_id, uint> getRaftSpecCountMap(
        solver::raft_global_id ridx,
        solver::spec_global_id spec_gidx) const;

    uint getRaftSpecCount(solver::raft_global_id ridx, solver::spec_global_id spec_gidx) const;

    const std::vector<triangle_global_id>& getRaftOverlap(solver::raft_global_id raft_gidx,
                                                          triangle_global_id tri_gidx) const {
        // THis is quite well protected, just throwing an error if
        // keys don't exist, so no more error checking really necessary
        return pTriOverlap.at(raft_gidx).at(tri_gidx);
    }

    void addTriRaftsRefs(triangle_global_id centraltri_gidx, Raft* raft);
    void removeTriRaftRefs(triangle_global_id centraltri_gidx, Raft* raft);

    // Tries to find a position for a raft. Returns false if can't be found (after
    // 10 attempts)
    bool getPos(math::position_abs* rpos,
                solver::Raftdef* raftdef,
                triangle_global_id central_tri_gidx,
                Raft* raft = nullptr) const;

    // Do Raft diffusion.
    void runRaft(double dt);

    inline TetVesicleVesRaft* solverVesRaft() const noexcept {
        return pVesRaft;
    }
    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    solver::Patchdef* pPatchdef;
    double pArea;

    tetmesh::Tetmesh* pMesh;

    // Just easier to store the solver for all the raft KProc stuff
    TetVesicleVesRaft* pVesRaft{nullptr};

    // Store a pointer to the RNG for convenience
    rng::RNGptr pRNG;

    util::strongid_vector<triangle_local_id, TriVesRaft*> pTris;

    std::map<triangle_global_id, triangle_local_id> pTriidcs_G_to_L;
    std::map<triangle_local_id, triangle_global_id> pTriidcs_L_to_G;

    // The rafts are stored as a map of raft type by index to a
    // list of rafts of that type. List should be more
    // efficient than a vector for insertions and removals
    std::map<solver::raft_global_id, std::list<Raft*>> pRafts;

    // Raft overlap structure raft_indx > tri_central_index > overlap_tris
    // That is for any raft index (key to map) gives a vector length of all tris,
    // each element giving a list of POSSIBLE overlap whenever raft is in that
    // triangle
    std::map<solver::raft_global_id, std::map<triangle_global_id, std::vector<triangle_global_id>>>
        pTriOverlap;

    // The zones, with list of triangles belonging to that zone
    std::map<std::string, std::vector<TriVesRaft*>> pZones;

    // Endocytosis index to map of zone and local (TetVesicleRDEF) endocytosis
    // object added to that zone
    std::map<solver::endocytosis_local_id, std::map<std::string, Endocytosis*>> pEndosMap;

    // Convenient also to store a simple vector
    std::vector<Endocytosis*> pEndosVec;
};

}  // namespace steps::mpi::tetvesicle
