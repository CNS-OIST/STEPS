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
#include "mpi/tetvesicle/kproc.hpp"
#include "mpi/tetvesicle/raft.hpp"
#include "mpi/tetvesicle/raftproxy.hpp"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"

namespace steps::mpi::tetvesicle {

// Forward declarations
class TetRDEF;
class TriRDEF;
class SReac;
class Exocytosis;
class SDiff;
class TetVesicleRDEF;
class VDepSReac;
class GHKcurr;
class RaftGen;
class PatchRDEF;

////////////////////////////////////////////////////////////////////////////////

// Auxiliary declarations.
typedef TriRDEF* TriRDEFP;
typedef std::vector<TriRDEF*> TriRDEFPVec;
typedef TriRDEFPVec::iterator TriRDEFPVecI;
typedef TriRDEFPVec::const_iterator TriRDEFPVecCI;

////////////////////////////////////////////////////////////////////////////////

class TriRDEF {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    TriRDEF(triangle_global_id idx,
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
            const math::point3d& trinorm,
            int rank,
            int host_rank);
    ~TriRDEF();

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
    void setInnerTet(TetRDEF* t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(TetRDEF* t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, TriRDEF* t);

    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(TetVesicleRDEF* tex, bool efield = false);

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

    inline void setPatchRDEF(PatchRDEF* patch) {
        AssertLog(pPatchRDEF == nullptr);
        pPatchRDEF = patch;
    }

    inline PatchRDEF* patchRDEF() const {
        AssertLog(pPatchRDEF != nullptr);
        return pPatchRDEF;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline double area() const noexcept {
        return pArea;
    }

    inline TetRDEF* iTet() const noexcept {
        return pInnerTet;
    }

    inline TetRDEF* oTet() const noexcept {
        return pOuterTet;
    }

    inline TriRDEF* nextTri(uint i) const {
        AssertLog(i < 3);
        return pNextTri[i];
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

    // Set whether a direction is a diffusion boundary
    void setSDiffBndDirection(uint i);

    inline bool getSDiffBndDirection(uint idx) const noexcept {
        return pSDiffBndDirection[idx];
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EFIELD
    ////////////////////////////////////////////////////////////////////////

    // Local index of GHK current given
    void incECharge(solver::ghkcurr_local_id lidx, int charge);

    // Should be called at the beginning of every EField time-step
    void resetECharge(double dt, double efdt);

    // reset the Ohmic current opening time integral info, also should be
    // called just before commencing or just after completing an EField dt
    void resetOCintegrals();

    double computeI(double v, double dt, double simtime, double efdt);

    double getOhmicI(double v, double dt) const;
    double getOhmicI(solver::ohmiccurr_local_id lidx, double v, double dt) const;

    double getGHKI(double dt) const;
    double getGHKI(solver::ghkcurr_local_id lidx, double dt) const;

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }

    void setCount(solver::spec_local_id lidx, uint count, double period = 0.0);
    void incCount(solver::spec_local_id lidx,
                  int inc,
                  double period = 0.0,
                  bool local_change = false);

    static const uint CLAMPED = 1;

    inline bool clamped(solver::spec_local_id lidx) const noexcept {
        return pPoolFlags[lidx] & CLAMPED;
    }
    void setClamped(solver::spec_local_id lidx, bool clamp);

    // Set a channel state relating to an ohmic current change.
    // 0th argument is oc local index, 1st argument is the local index
    // of the related channel state
    void setOCchange(solver::ohmiccurr_local_id oclidx,
                     solver::spec_local_id slidx,
                     double dt,
                     double simtime);

    // Set/get the reversal potential of an ohmic current
    void setOCerev(solver::ohmiccurr_local_id oclidx, double erev);
    double getOCerev(solver::ohmiccurr_local_id oclidx) const;

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<KProc*>::const_iterator kprocBegin() const noexcept {
        return pKProcs.begin();
    }
    inline std::vector<KProc*>::const_iterator kprocEnd() const noexcept {
        return pKProcs.end();
    }
    inline std::vector<KProc*> const& kprocs() const noexcept {
        return pKProcs;
    }
    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }

    inline KProc* getKProc(uint lidx) const {
        if (hostRank != myRank)
            return nullptr;
        AssertLog(lidx < pKProcs.size());
        return pKProcs[lidx];
    }

    SReac& sreac(solver::sreac_local_id lidx) const;
    RaftGen& raftgen(solver::raftgen_local_id lidx) const;
    SDiff& sdiff(solver::surfdiff_local_id lidx) const;
    VDepSReac& vdepsreac(solver::vdepsreac_local_id lidx) const;
    GHKcurr& ghkcurr(solver::ghkcurr_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////

    inline math::point3d const& position() const noexcept {
        return pPosition;
    }

    inline math::point3d const& norm() const noexcept {
        return pNorm;
    }

    ////////////////////////////////////////////////////////////////////////
    // RAFT-RELATED
    ////////////////////////////////////////////////////////////////////////

    void createRaftProxyref(solver::Raftdef* raftdef, solver::raft_individual_id raft_unique_id);

    inline std::map<solver::raft_individual_id, RaftProxy*> const& getRaftProxyrefs()
        const noexcept {
        return pRaftProxyrefs;
    }

    RaftProxy* getRaftProxyref(solver::raft_individual_id raft_uidx);

    void clearRaftProxyrefs();

    ///////////////////// ADDED FOR MPI STEPS ////////////////

    // setup dependence for KProcs in this subvolume
    void setupDeps();

    inline uint getStartKProcIdx() const noexcept {
        return startKProcIdx;
    }

    // check if kp_lidx in this vol depends on spec_gidx in WMVol kp_container
    bool KProcDepSpecTet(uint kp_lidx, TetRDEF* kp_container, solver::spec_global_id spec_gidx);
    // check if kp_lidx in this vol depends on spec_gidx in Tri kp_container
    bool KProcDepSpecTri(uint kp_lidx, TriRDEF* kp_container, solver::spec_global_id spec_gidx);
    // check if kp_lidx in this vol depends on spec_gidx in Tri kp_container
    // specifically on raft surface (any??)
    bool KProcDepSpecTriRaftSurface(uint kp_lidx,
                                    TriRDEF* kp_container,
                                    solver::spec_global_id spec_gidx);

    bool getInHost() const;
    void setHost(int host, int rank);
    int getHost() const noexcept {
        return hostRank;
    }

    void setSolverRDEF(TetVesicleRDEF* solver);
    TetVesicleRDEF* solverRDEF() const;

    double getPoolOccupancy(solver::spec_local_id lidx);
    double getLastUpdate(solver::spec_local_id lidx);
    void resetPoolOccupancy();

    std::vector<KProc*> const& getSpecUpdKProcs(solver::spec_local_id slidx) const noexcept {
        return localSpecUpdKProcs[slidx];
    }

    void setupBufferLocations();

    ////////////////////////////////////////////////////////////////////////

    // Flag a raftgen to be applied at next opportunity. Remove the species from
    // the triangle so they can't get used by anything else.
    void applyRaftGen(solver::RaftGendef* raftgen, double period);

    inline const std::map<solver::raftgen_global_id, uint>& getAppliedRaftGens() const noexcept {
        return pAppliedRaftgens;
    }

    inline void clearAppliedRaftGens() {
        pAppliedRaftgens.clear();
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    triangle_global_id pIdx;

    solver::Patchdef* pPatchdef;

    double pArea;

    /// Pointers to neighbouring tetrahedra.
    TetRDEF* pInnerTet;
    TetRDEF* pOuterTet;

    // Indices of two neighbouring tets; UNKNOWN_TET if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    tetrahedron_global_id pTets[2];

    // Indices of neighbouring triangles.
    triangle_global_id pTris[3];

    /// POinters to neighbouring triangles
    TriRDEF* pNextTri[3];

    bool hasEfield{};

    // Neighbour information- needed for surface diffusion
    double pLengths[3];
    double pDist[3];

    bool pSDiffBndDirection[3];

    /// Numbers of molecules -- stored as machine word integers.
    util::strongid_vector<solver::spec_local_id, uint> pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    util::strongid_vector<solver::spec_local_id, uint> pPoolFlags;

    /// The kinetic processes.
    std::vector<KProc*> pKProcs;
    uint startKProcIdx{};
    uint nKProcs{};

    /// For the EFIELD calculation. An integer storing the amount of
    /// elementary charge from inner tet to outer tet (positive if
    /// net flux is positive, negative if net flux is negative) for
    /// one EField time-step.
    // NOTE: Now arrays so as to separate into different GHK currs,
    // for data access
    util::strongid_vector<solver::ghkcurr_local_id, int> pECharge;

    // to store the latest ECharge, so that the info is available to solver
    util::strongid_vector<solver::ghkcurr_local_id, int> pECharge_last;
    util::strongid_vector<solver::ghkcurr_local_id, int> pECharge_accum;
    double pECharge_last_dt;
    double pECharge_accum_dt;

    // Store the Ohmic currents' channel's opening information by OC local indices
    // and the time since the related channel state changed it's number
    // The pOCchan_timeintg stores number of channel open * opening time
    // so at the end of the step this number/Efield dt will give the
    // mean number of channels open
    util::strongid_vector<solver::ohmiccurr_local_id, double> pOCchan_timeintg;

    util::strongid_vector<solver::ohmiccurr_local_id, double> pOCtime_upd;

    // Store reversal potential here to enable modification within API
    std::map<solver::ohmiccurr_local_id, double> pERev;

    math::point3d pPosition;

    math::point3d pNorm;

    // Now store a set of overlap rafts- useful information for updates
    std::map<solver::raft_individual_id, RaftProxy*> pRaftProxyrefs;

    PatchRDEF* pPatchRDEF;

    /////////////// MPI FUNCTIONALITY
    int hostRank;
    int myRank;
    TetVesicleRDEF* pRDEF{nullptr};

    // there is no need for sync tri, because no kproc in other surface or volume
    // depends on molecule changes of this surface
    std::set<int> syncHosts;

    util::strongid_vector<solver::spec_local_id, double> pPoolOccupancy;
    /// Structure to store time since last update, used to calculate occupancy
    util::strongid_vector<solver::spec_local_id, double> pLastUpdate;

    util::strongid_vector<solver::spec_local_id, uint> bufferLocations;
    util::strongid_vector<solver::spec_local_id, std::vector<KProc*>> localSpecUpdKProcs;

    ////////////////////////////////////////////////////////////////////////

    std::map<solver::raftgen_global_id, uint> pAppliedRaftgens;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::mpi::tetvesicle
