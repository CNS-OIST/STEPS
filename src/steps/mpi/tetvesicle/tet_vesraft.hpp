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
#include "math/point.hpp"
#include "solver/compdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"
#include "util/error.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations
class TetVesRaft;
class TriVesRaft;
class TetVesicleVesRaft;
class Vesicle;
class VesBind;
class CompVesRaft;

// Auxiliary declarations.
// typedef TetVesRaft *TetVesRaftP;
typedef std::vector<TetVesRaft*> TetVesRaftPVec;
typedef TetVesRaftPVec::iterator TetVesRaftPVecI;
typedef TetVesRaftPVec::const_iterator TetVesRaftPVecCI;

////////////////////////////////////////////////////////////////////////////////

class TetVesRaft {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    TetVesRaft(tetrahedron_global_id idx,
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
               math::point3d baryc);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////

    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////

    void reset();

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////

    inline solver::Compdef* compdef() const noexcept {
        return pCompdef;
    }

    inline tetrahedron_global_id idx() const noexcept {
        return pIdx;
    }

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set pointer to the next neighbouring tetrahedron.
    ///
    void setNextTet(uint i, TetVesRaft* t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, TriVesRaft* t);

    inline void setCompVesRaft(CompVesRaft* comp) {
        AssertLog(pCompVesRaft == nullptr);
        pCompVesRaft = comp;
    }

    inline CompVesRaft* getCompVesRaft() const {
        AssertLog(pCompVesRaft != nullptr);
        return pCompVesRaft;
    }

    ////////////////////////////////////////////////////////////////////////
    // VESICLE-RELATED
    ////////////////////////////////////////////////////////////////////////

    // Add a vesicle reference (vesicle overlaps this tet)
    inline void addVesref(Vesicle* ves) noexcept {
        pVesrefs.insert(ves);
    }

    // Remove a vesicle reference (vesicle has left this tet)
    inline void removeVesref(Vesicle* ves) noexcept {
        pVesrefs.erase(ves);
    }

    inline std::set<Vesicle*, util::DerefPtrLess<Vesicle>> const& getVesrefs() const noexcept {
        return pVesrefs;
    }

    // can be used to either increase or reduce
    void changeOverlap(double overlap);

    inline bool isFullOverlap() const noexcept {
        return pOverlap == pVol;
    }

    inline double getOverlap() const noexcept {
        return pOverlap;
    }

    ////////////////////////////////////////////////////////////////////////
    // SPECIES-RELATED
    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    void setCount(solver::spec_local_id lidx, uint count);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume without considering vesicle overlap
    inline double staticVol() const noexcept {
        return pVol;
    }

    /// Get the volume. (Effecive) volume changes with vesicle occupancy.
    ///
    inline double vol() const noexcept {
        return pReducedVol;
    }

    // Return the barycentre
    inline math::point3d& position() {
        return pPosition;
    }

    ///////////////////// Neighbouring tetrahedrons ////////////////////////

    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline TetVesRaft* nextTet(uint i) const noexcept {
        return pNextTet[i];
    }

    inline tetrahedron_global_id tet(uint t) const noexcept {
        return pTets[t];
    }

    /// Find the direction index towards a neighbor tetrahedron.
    ///
    int getTetDirection(tetrahedron_global_id tidx);

    /// Get the distance to the centroid of the next neighbouring
    /// tetrahedron.
    ///
    inline double dist(uint i) const noexcept {
        return pDist[i];
    }
    ///////////////////// Neighbouring triangles ////////////////////////

    /// Get pointer to the next neighbouring triangle.
    ///
    inline TriVesRaft* nextTri(uint i) const {
        AssertLog(i < 4);
        return pNextTris[i];
    }

    inline auto nexttriBegin() const noexcept {
        return pNextTris.begin();
    }
    inline auto nexttriEnd() const noexcept {
        return pNextTris.end();
    }
    inline const auto& nexttris() const noexcept {
        return pNextTris;
    }

    /// Get the area of a boundary triangle.
    ///
    inline double area(uint i) const noexcept {
        return pAreas[i];
    }

    //////////////////////////// MPI STEPS /////////////////////////////////

    void setSolverVesRaft(TetVesicleVesRaft* solver) noexcept {
        pVesRaft = solver;
    }

    TetVesicleVesRaft* solverVesRaft() const noexcept {
        return pVesRaft;
    };

    ////////////////////////////////////////////////////////////////////////

  private:
    //////////////////// This stuff copied from WmVol /////////////////

    tetrahedron_global_id pIdx;

    solver::Compdef* pCompdef;

    double pVol;

    // Overlap with vesicles
    double pOverlap;

    double pReducedVol;

    /// Numbers of molecules -- stored as uint.
    util::strongid_vector<solver::spec_local_id, uint> pPoolCount;

    math::point3d pPosition;

    CompVesRaft* pCompVesRaft{};

    std::set<Vesicle*, util::DerefPtrLess<Vesicle>> pVesrefs;

    ////////////////////////////////////////////////////////////////////////

    // Indices of neighbouring tetrahedra.
    std::array<tetrahedron_global_id, 4> pTets;

    /// Pointers to neighbouring tetrahedra.
    std::array<TetVesRaft*, 4> pNextTet;

    // The connected patch triangles.
    std::array<TriVesRaft*, 4> pNextTris{};

    std::array<double, 4> pAreas;
    std::array<double, 4> pDist;

    ///////// MPI STEPS ////////////////////////////////////////////////////

    TetVesicleVesRaft* pVesRaft{nullptr};

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
