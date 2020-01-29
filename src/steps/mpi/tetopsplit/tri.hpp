/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


#ifndef STEPS_MPI_TETOPSPLIT_TRI_HPP
#define STEPS_MPI_TETOPSPLIT_TRI_HPP 1

// STL headers.
#include <cassert>
#include <vector>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/solver/types.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations

class Tet;
class WmVol;
class Tri;
class SReac;
class SDiff;
class TetOpSplitP;
class VDepTrans;
class VDepSReac;
class GHKcurr;

////////////////////////////////////////////////////////////////////////////////

// Auxiliary declarations.
typedef Tri *                           TriP;
typedef std::vector<TriP>               TriPVec;
typedef TriPVec::iterator               TriPVecI;
typedef TriPVec::const_iterator         TriPVecCI;

////////////////////////////////////////////////////////////////////////////////

class Tri
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Tri(triangle_id_t idx, steps::solver::Patchdef *patchdef, double area,
        double l0, double l1, double l2, double d0, double d1, double d2,
        tetrahedron_id_t tetinner, tetrahedron_id_t tetouter,
        triangle_id_t tri0, triangle_id_t tri1, triangle_id_t tri2, int rank, int host_rank);
    virtual ~Tri();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set pointer to the 'inside' neighbouring tetrahedron.
    ///
    void setInnerTet(smtos::WmVol * t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(smtos::WmVol * t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, smtos::Tri * t);
    
    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(smtos::TetOpSplitP * tex, bool efield = false);

    /// Set all pool flags and molecular populations to zero.
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Patchdef * patchdef() const noexcept
    { return pPatchdef; }

    inline triangle_id_t idx() const noexcept
    { return pIdx; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline double area() const noexcept
    { return pArea; }

    inline smtos::WmVol * iTet() const noexcept
    { return pInnerTet; }

    inline smtos::WmVol * oTet() const noexcept
    { return pOuterTet; }

    inline smtos::Tri * nextTri(uint i) const
    {
        AssertLog(i < 3);
        return pNextTri[i];
    }

    inline triangle_id_t tri(uint t) const noexcept
    { return pTris[t]; }

    /// Get the length of a boundary bar.
    ///
    inline double length(uint i) const noexcept
    { return pLengths[i]; }

    /// Get the distance to the centroid of the next neighbouring
    /// triangle.
    ///
    inline double dist(uint i) const noexcept
    { return pDist[i]; }

    inline tetrahedron_id_t tet(uint t) const noexcept
    { return pTets[t]; }
    
    /// Find the direction index towards a neighbor triangle.
    ///
    int getTriDirection(triangle_id_t tidx);

    ////////////////////////////////////////////////////////////////////////

    // Set whether a direction is a surface diffusion boundary
    void setSDiffBndDirection(uint i);

    inline bool getSDiffBndDirection(uint idx) const noexcept
    { return pSDiffBndDirection[idx]; }

    /////////////////////////// Dependency ////////////////////////////////
    inline uint getStartKProcIdx() const noexcept
    {return startKProcIdx;}
    
    // setup dependence for KProcs in this subvolume
    void setupDeps();
    
    // check if kp_lidx in this vol depends on spec_gidx in WMVol kp_container
    virtual bool KProcDepSpecTet(uint kp_lidx, WmVol* kp_container, uint spec_gidx);
    // check if kp_lidx in this vol depends on spec_gidx in Tri kp_container
    virtual bool KProcDepSpecTri(uint kp_lidx, Tri* kp_container, uint spec_gidx);
    
    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EFIELD
    ////////////////////////////////////////////////////////////////////////

    // Local index of GHK current given
    void incECharge(uint lidx, int charge);

    // Should be called at the beginning of every EField time-step
    void resetECharge();

    // reset the Ohmic current opening time integral info, also should be
    // called just before commencing or just after completing an EField dt
    void resetOCintegrals();

    double computeI(double v, double dt, double simtime);

    double getOhmicI(double v, double dt) const;
    double getOhmicI(uint lidx, double v,double dt) const;

    double getGHKI( double dt) const;
    double getGHKI(uint lidx, double dt) const;

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    inline uint * pools() const
    { return pPoolCount; }
    void setCount(uint lidx, uint count, double period = 0.0);
    void incCount(uint lidx, int inc, double period = 0.0, bool local_change = false);


    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const noexcept
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    // Set a channel state relating to an ohmic current change.
    // 0th argument is oc local index, 1st argument is the local index
    // of the related channel state
    void setOCchange(uint oclidx, uint slidx, double dt, double simtime);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<smtos::KProc *>::const_iterator kprocBegin() const noexcept
    { return pKProcs.begin(); }
    inline std::vector<smtos::KProc *>::const_iterator kprocEnd() const noexcept
    { return pKProcs.end(); }
    inline std::vector<smtos::KProc *> & kprocs() noexcept
    { return pKProcs; }
    inline uint countKProcs() const noexcept
    { return nKProcs; }

    inline KProc * getKProc(uint lidx) const
    {
        if (hostRank != myRank) return nullptr;
        AssertLog(lidx < pKProcs.size());
        return pKProcs[lidx];
    }
    
    smtos::SReac * sreac(uint lidx) const;
    smtos::SDiff * sdiff(uint lidx) const;
    smtos::VDepTrans * vdeptrans(uint lidx) const;
    smtos::VDepSReac * vdepsreac(uint lidx) const;
    smtos::GHKcurr * ghkcurr(uint lidx) const;

    ////////////////////////////////////////////////////////////////////////

    bool getInHost();
    void setHost(int host, int rank);
    int getHost() {return hostRank;}
    void setSolver(steps::mpi::tetopsplit::TetOpSplitP* sol);
    steps::mpi::tetopsplit::TetOpSplitP* solver();
    
    //void sendSyncPools();
    //void recvSyncPools(int source);
	double getPoolOccupancy(uint lidx);
	double getLastUpdate(uint lidx);
	void resetPoolOccupancy();

    std::vector<smtos::KProc*> const & getSpecUpdKProcs(uint slidx);
    
    void repartition(smtos::TetOpSplitP * tex, int rank, int host_rank);
    void setupBufferLocations();
private:

    ////////////////////////////////////////////////////////////////////////

    triangle_id_t                       pIdx;

    steps::solver::Patchdef           * pPatchdef;

    double                              pArea;

    double                              pLengths[3];
    double                              pDist[3];

  /// Pointers to neighbouring tetrahedra.
    smtos::WmVol                       * pInnerTet{nullptr};
    smtos::WmVol                       * pOuterTet{nullptr};

    // Indices of two neighbouring tets; -1 if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    tetrahedron_id_t                     pTets[2];

    // Indices of neighbouring triangles.
    triangle_id_t                        pTris[3];

    /// POinters to neighbouring triangles
    smtos::Tri                         * pNextTri[3];

    bool                                hasEfield;

    // Neighbour information- needed for surface diffusion

    bool                                pSDiffBndDirection[3];

    /// Numbers of molecules -- stored as machine word integers.
    uint                              * pPoolCount{nullptr};
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags{nullptr};

    /// The kinetic processes.
    std::vector<smtos::KProc *>          pKProcs;
    uint                                 startKProcIdx;
    uint                                 nKProcs;
    /// For the EFIELD calculation. An integer storing the amount of
    /// elementary charge from inner tet to outer tet (positive if
    /// net flux is positive, negative if net flux is negative) for
    /// one EField time-step.
    // NOTE: Now arrays so as to separate into different GHK currs,
    // for data access
    int                                  * pECharge{nullptr};

    // to store the latest ECharge, so that the info is available to solver
    int                               * pECharge_last{nullptr};


    // Store the Ohmic currents' channel's opening information by OC local indices
    // and the time since the related channel state changed it's number
    // The pOCchan_timeintg stores number of channel open * opening time
    // so at the end of the step this number/Efield dt will give the
    // mean number of channels open
    double                               * pOCchan_timeintg{nullptr};

    double                               * pOCtime_upd{nullptr};
    ///////////////MPI STUFFS
    int                                     hostRank;
    int                                     myRank;
    steps::mpi::tetopsplit::TetOpSplitP   * pSol;
    
	// there is no need for sync tri, because no kproc in other surface or volume depends on molecule changes
	// of this surface
    std::set<int>                           syncHosts;

    double                            * pPoolOccupancy{nullptr};
    /// Structure to store time since last update, used to calculate occupancy
    double 							  *	pLastUpdate{nullptr};
    
    std::vector<uint>                   bufferLocations;
    std::vector<std::vector<smtos::KProc *>> localSpecUpdKProcs;
    

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_TRI_HPP

// END
