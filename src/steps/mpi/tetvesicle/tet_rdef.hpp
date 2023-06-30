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
#include "mpi/tetvesicle/kproc.hpp"
#include "solver/compdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations
class TetRDEF;
class Diff;
class TriRDEF;
class Reac;
class TetVesicleRDEF;
class VesProxy;
class VesBind;
class CompRDEF;

// Auxiliary declarations.
typedef TetRDEF* TetRDEFP;
typedef std::vector<TetRDEF*> TetRDEFPVec;
typedef TetRDEFPVec::iterator TetRDEFPVecI;
typedef TetRDEFPVec::const_iterator TetRDEFPVecCI;

////////////////////////////////////////////////////////////////////////////////

// Used to store information about new link specs, to be aded to vesicles
// on next update

struct LinkSpec_newpair {
    solver::vesicle_individual_id vesicle1_individual_id;
    solver::vesicle_individual_id vesicle2_individual_id;

    solver::linkspec_global_id linkspec1_global_id;
    solver::linkspec_global_id linkspec2_global_id;

    solver::linkspec_individual_id linkspec1_individual_id;
    solver::linkspec_individual_id linkspec2_individual_id;

    math::position_abs linkspec1_pos_absolute;
    math::position_abs linkspec2_pos_absolute;

    double min_length;
    double max_length;
};

////////////////////////////////////////////////////////////////////////////////

class TetRDEF {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    TetRDEF(tetrahedron_global_id idx,
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
            const math::point3d& baryc,
            int rank,
            int host_rank);
    ~TetRDEF();

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
    void setNextTet(uint i, TetRDEF* t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, TriRDEF* t);

    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(TetVesicleRDEF* tex);

    inline void setCompRDEF(CompRDEF* comp) {
        AssertLog(pCompRDEF == nullptr);
        pCompRDEF = comp;
    }

    inline CompRDEF* getCompRDEF() const {
        AssertLog(pCompRDEF != nullptr);
        return pCompRDEF;
    }

    ////////////////////////////////////////////////////////////////////////
    // VESICLE-RELATED
    ////////////////////////////////////////////////////////////////////////

    // Add a vesicle reference (vesicle overlaps this tet)
    void createVesProxyref(solver::Vesicledef* vesdef,
                           solver::vesicle_individual_id ves_unique_id,
                           math::position_abs ves_pos,
                           bool contains_link);

    VesProxy* getVesProxyref(solver::vesicle_individual_id veuidx);

    inline std::map<solver::vesicle_individual_id, VesProxy*> const& getVesProxyrefs()
        const noexcept {
        return pVesProxyrefs;
    }

    void clearVesProxyrefs();

    // Set the overlap volume with vesicles
    void setOverlap(double overlap);

    ////////////////////////////////////////////////////////////////////////
    // SPECIES-RELATED
    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    void setCount(solver::spec_local_id lidx, uint count, double period = 0.0);
    void incCount(solver::spec_local_id lidx,
                  int inc,
                  double period = 0.0,
                  bool local_change = false);

    // The concentration of species global index gidx in MOL PER l
    double conc(solver::spec_global_id gidx) const;

    static const uint CLAMPED = 1;

    inline bool clamped(solver::spec_local_id lidx) const noexcept {
        return (pPoolFlags[lidx] & CLAMPED) != 0;
    }

    void setClamped(solver::spec_local_id lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////
    // KPROCS
    ////////////////////////////////////////////////////////////////////////

    inline std::vector<KProc*>::const_iterator kprocBegin() const noexcept {
        return pKProcs.begin();
    }
    inline std::vector<KProc*>::const_iterator kprocEnd() const noexcept {
        return pKProcs.end();
    }
    inline uint getStartKProcIdx() const noexcept {
        return startKProcIdx;
    }
    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }
    inline std::vector<KProc*>& kprocs() noexcept {
        return pKProcs;
    }

    // Added for MPI
    inline KProc* getKProc(uint lidx) {
        if (hostRank != myRank) {
            return nullptr;
        }
        AssertLog(lidx < pKProcs.size());
        return pKProcs[lidx];
    }

    Reac& reac(solver::reac_local_id lidx) const;

    Diff& diff(solver::diff_local_id lidx) const;

    VesBind& vesbind(solver::vesbind_local_id lidx) const;

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
        double vol = pVol - pOverlap;
        // Never allow zero volume since that can give NANs in diff.
        // THis is an issue if vesicles cover tets before they have
        // evacuated which can happen in very fine meshes or bad models
        if (vol <= 0.0) {
            return solver::TINY_VOLUME;
        } else {
            return vol;
        }
    }

    // Return the barycentre
    // inline math::point3d &position() { return pPosition; }

    ///////////////////// Neighbouring tetrahedrons ////////////////////////

    /// Get pointer to the next neighbouring tetrahedron.
    ///
    inline TetRDEF* nextTet(uint i) const noexcept {
        return pNextTet[i];
    }

    const auto& nextTets() const noexcept {
        return pNextTet;
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
    inline TriRDEF* nextTri(uint i) const {
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

    ///////////////////////// Diffusion Boundary //////////////////////////

    // Set whether a direction is a diffusion boundary
    void setDiffBndDirection(uint i);

    inline bool getDiffBndDirection(uint idx) const noexcept {
        return pDiffBndDirection[idx];
    }

    //////////////////////////// MPI STEPS /////////////////////////////////

    bool getInHost() const;
    inline int getHost() const {
        return hostRank;
    }
    void setHost(int host, int rank);

    void setSolverRDEF(TetVesicleRDEF* solver);
    TetVesicleRDEF* solverRDEF() const;

    void setupBufferLocations();

    // setup dependence for KProcs in this subvolume
    void setupDeps();

    // check if kp_lidx in this vol depends on spec_gidx in Tet kp_container
    bool KProcDepSpecTet(uint kp_lidx, TetRDEF* kp_container, solver::spec_global_id spec_gidx);

    // check if kp_lidx in this vol depends on linkspec_gidx in Tet kp_container
    // (linkspecs always on ves surface)
    bool KProcDepLinkSpecTetVesSurface(uint kp_lidx,
                                       TetRDEF* kp_container,
                                       solver::linkspec_global_id linkspec_gidx);

    // An extension to KProcDepSpecTet to check if kp_lidx in this vol depends
    // on spec_gidx in Tet kp_container but specifically on (any??) vesicle
    // surface
    bool KProcDepSpecTetVesSurface(uint kp_lidx,
                                   TetRDEF* kp_container,
                                   solver::spec_global_id spec_gidx);

    // check if kp_lidx in this vol depends on spec_gidx in Tri kp_container
    bool KProcDepSpecTri(uint kp_lidx,
                         TriRDEF* kp_container,
                         solver::spec_global_id spec_gidx) const;

    inline double getPoolOccupancy(solver::spec_local_id lidx) const {
        AssertLog(lidx < compdef()->countSpecs());
        return pPoolOccupancy[lidx];
    }

    double getLastUpdate(solver::spec_local_id lidx) const {
        AssertLog(lidx < compdef()->countSpecs());
        return pLastUpdate[lidx];
    }

    void resetPoolOccupancy();

    std::vector<KProc*> const& getSpecUpdKProcs(solver::spec_local_id slidx) const noexcept {
        return localSpecUpdKProcs[slidx.get()];
    }


    solver::linkspec_individual_id getNextLinkSpecUniqueIndex() const;

    ////////////////////////////////////////////////////////////////////////

    void addNewLinkedSpecs(solver::linkspec_global_id linkspec1_global_id,
                           solver::linkspec_global_id linkspec2_global_id,
                           solver::linkspec_individual_id linkspec1,
                           solver::linkspec_individual_id linkspec2,
                           VesProxy* ves1,
                           VesProxy* ves2,
                           math::position_abs linkspec1_pos_abs,
                           math::position_abs linkspec2_pos_abs,
                           double max_length,
                           double min_length);

    inline bool haveNewLinkedSpecs() {
        return !pLinkedSpecsUpd.empty();
    }

    inline std::vector<LinkSpec_newpair> const& getNewLinkedSpecs() {
        return pLinkedSpecsUpd;
    }

    inline void clearNewLinkedSpecs() {
        pLinkedSpecsUpd.clear();
    }

    ////////////////////////////////////////////////////////////////////////

    // protected:
  private:
    /// The kinetic processes.
    std::vector<KProc*> pKProcs;
    uint startKProcIdx{};
    uint nKProcs;

    //////////////////// This stuff copied from WmVol /////////////////

    tetrahedron_global_id pIdx;

    solver::Compdef* pCompdef;

    double pVol;

    // Overlap with vesicles
    double pOverlap;

    /// Numbers of molecules -- stored as uint.
    steps::util::strongid_vector<solver::spec_local_id, uint> pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    steps::util::strongid_vector<solver::spec_local_id, uint> pPoolFlags;

    CompRDEF* pCompRDEF;

    // Now store a set of overlap vesicles- useful information for updates
    std::map<solver::vesicle_individual_id, VesProxy*> pVesProxyrefs;

    // Have to store any newly created linked species to tell comp at the next
    // opportunity
    std::vector<LinkSpec_newpair> pLinkedSpecsUpd;  // the pair

    ////////////////////////////////////////////////////////////////////////

    // Indices of neighbouring tetrahedra.
    tetrahedron_global_id pTets[4];

    /// Pointers to neighbouring tetrahedra.
    std::array<TetRDEF*, 4> pNextTet{};

    // The connected patch triangles.
    std::array<TriRDEF*, 4> pNextTris{};

    std::array<double, 4> pAreas;
    std::array<double, 4> pDist;

    std::array<bool, 4> pDiffBndDirection{};

    ///////// MPI STEPS ////////////////////////////////////////////////////
    int myRank;
    int hostRank;
    TetVesicleRDEF* pRDEF{nullptr};

    /// Structure to store occupancy over the update period
    steps::util::strongid_vector<solver::spec_local_id, double> pPoolOccupancy;
    /// Structure to store time since last update, used to calculate occupancy
    steps::util::strongid_vector<solver::spec_local_id, double> pLastUpdate;

    /// location of where the change of this species is stored in  the solver
    /// buffer
    util::strongid_vector<solver::spec_local_id, uint> bufferLocations;
    // local kprocs update list when spec is changed
    std::vector<std::vector<KProc*>> localSpecUpdKProcs;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
