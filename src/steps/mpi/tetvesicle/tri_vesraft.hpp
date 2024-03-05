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
#include <vector>

// STEPS headers.
#include "mpi/tetvesicle/raft.hpp"
#include "solver/patchdef.hpp"
#include "util/error.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations

class TetVesRaft;
class TriVesRaft;
class TetVesicleVesRaft;

////////////////////////////////////////////////////////////////////////////////

class TriVesRaft {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    TriVesRaft(triangle_global_id idx,
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
               const math::point3d& trinorm);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set pointer to the 'inside' neighbouring tetrahedron.
    ///
    void setInnerTet(TetVesRaft* t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(TetVesRaft* t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, TriVesRaft* t);

    /// Set all pool flags and molecular populations to zero.
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////

    inline solver::Patchdef* patchdef() const noexcept {
        return pPatchdef;
    }

    inline triangle_global_id idx() const noexcept {
        return pIdx;
    }

    inline void setPatchVesRaft(PatchVesRaft* patch) {
        AssertLog(pPatchVesRaft == nullptr);
        pPatchVesRaft = patch;
    }

    inline PatchVesRaft* patchVesRaft() const {
        AssertLog(pPatchVesRaft != nullptr);
        return pPatchVesRaft;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline double area() const noexcept {
        return pArea;
    }

    inline TetVesRaft* iTet() const noexcept {
        return pInnerTet;
    }

    inline TetVesRaft* oTet() const noexcept {
        return pOuterTet;
    }

    inline TriVesRaft* nextTri(uint i) const {
        return pNextTri.at(i);
    }

    inline triangle_global_id tri(uint t) const noexcept {
        return pTris[t];
    }

    /// Get the length of a boundary bar.
    ///
    inline double length(uint i) const noexcept {
        return pLengths[i];
    }

    /// Get the distance to the centroid of the next neighbouring
    /// triangle.
    ///
    inline double dist(uint i) const noexcept {
        return pDist[i];
    }

    inline tetrahedron_global_id tet(uint t) const noexcept {
        return pTets[t];
    }

    /// Find the direction index towards a neighbor triangle.
    ///
    int getTriDirection(triangle_global_id tidx);

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    void setCount(solver::spec_local_id lidx, uint count);
    // void incCount(solver::spec_local_id lidx, int inc);

    ////////////////////////////////////////////////////////////////////////

    inline math::point3d const& position() noexcept {
        return pPosition;
    }

    inline math::point3d const& norm() noexcept {
        return pNorm;
    }

    ////////////////////////////////////////////////////////////////////////
    // RAFT-RELATED
    ////////////////////////////////////////////////////////////////////////

    inline void addRaftref(Raft* raft) noexcept {
        pRaftrefs.insert(raft);
    }
    inline void removeRaftref(Raft* raft) noexcept {
        pRaftrefs.erase(raft);
    }

    inline std::set<Raft*, util::DerefPtrLess<Raft>> const& getRaftrefs() const noexcept {
        return pRaftrefs;
    }

    inline const std::map<solver::raftgen_global_id, uint>& getAppliedRaftGens() const noexcept {
        return pAppliedRaftgens;
    }

    inline void clearAppliedRaftGens() {
        pAppliedRaftgens.clear();
    }

    inline void addRaftGen(solver::raftgen_global_id raftgengidx, uint count) {
        pAppliedRaftgens[raftgengidx] = count;
    }


    ///////////////////// ADDED FOR MPI STEPS ////////////////

    void setSolverVesRaft(TetVesicleVesRaft* solver);
    TetVesicleVesRaft* solverVesRaft() const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    triangle_global_id pIdx;

    solver::Patchdef* pPatchdef;

    double pArea;

    /// Pointers to neighbouring tetrahedra.
    TetVesRaft* pInnerTet;
    TetVesRaft* pOuterTet;

    // Indices of two neighbouring tets; UNKNOWN_TET if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    tetrahedron_global_id pTets[2];

    // Indices of neighbouring triangles.
    std::array<triangle_global_id, 3> pTris;

    /// POinters to neighbouring triangles
    std::array<TriVesRaft*, 3> pNextTri;

    // Neighbour information- needed for surface diffusion
    double pLengths[3];
    double pDist[3];

    /// Numbers of molecules -- stored as machine word integers.
    util::strongid_vector<solver::spec_local_id, uint> pPoolCount;

    math::point3d pPosition;

    math::point3d pNorm;

    // Now store a set of overlap rafts- useful information for updates
    std::set<Raft*, util::DerefPtrLess<Raft>> pRaftrefs;

    PatchVesRaft* pPatchVesRaft;

    std::map<solver::raftgen_global_id, uint> pAppliedRaftgens;

    /////////////// MPI FUNCTIONALITY
    TetVesicleVesRaft* pVesRaft{nullptr};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
