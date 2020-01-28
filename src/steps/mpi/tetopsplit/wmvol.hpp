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


#ifndef STEPS_MPI_TETOPSPLIT_WMVOL_HPP
#define STEPS_MPI_TETOPSPLIT_WMVOL_HPP 1

// STL headers.
#include <cassert>
#include <vector>
#include <set>
#include <fstream>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/compdef.hpp"
#include "steps/mpi/tetopsplit/kproc.hpp"
#include "steps/solver/types.hpp"
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace mpi {
 namespace tetopsplit {

////////////////////////////////////////////////////////////////////////////////

namespace smtos = steps::mpi::tetopsplit;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class WmVol;
class Tri;
class Reac;
class Tet;
class TetOpSplitP;

// Auxiliary declarations.
typedef WmVol *                           WmVolP;
typedef std::vector<WmVolP>               WmVolPVec;
typedef WmVolPVec::iterator               WmVolPVecI;
typedef WmVolPVec::const_iterator         WmVolPVecCI;

////////////////////////////////////////////////////////////////////////////////

// Base class for the tetrahedrons in the mesh. This allows for compartments to
// be described as a well-mixed volume or comprised of tetrahedrons in the
// reaction-diffusion solver. Of course, if a compartment is well-mixed,
// any diffusion rules are ignored.
//
class WmVol
{

public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    WmVol
      (
        tetrahedron_id_t idx, steps::solver::Compdef *cdef, double vol, int rank, int host_rank
      );

    virtual ~WmVol();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    virtual void checkpoint(std::fstream & cp_file);

    /// restore data
    virtual void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    virtual void setupKProcs(smtos::TetOpSplitP * tex);

    virtual void setNextTri(smtos::Tri *t);

    ////////////////////////////////////////////////////////////////////////

    virtual void reset();

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    inline steps::solver::Compdef * compdef() const noexcept
    { return pCompdef; }

    inline tetrahedron_id_t idx() const noexcept
    { return pIdx; }

    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume.
    ///
    inline double vol() const noexcept
    { return pVol; }


    ////////////////////////////////////////////////////////////////////////

    inline uint * pools() const
    { return pPoolCount; }
    virtual void setCount(uint lidx, uint count, double period = 0.0);
    virtual void incCount(uint lidx, int inc, double period = 0.0, bool local_change = false);

    // The concentration of species global index gidx in MOL PER l
    double conc(uint gidx) const;

    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const noexcept
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<smtos::Tri *>::const_iterator nexttriBegin() const noexcept
    { return pNextTris.begin(); }
    inline std::vector<smtos::Tri *>::const_iterator nexttriEnd() const noexcept
    { return pNextTris.end(); }
    inline const std::vector<smtos::Tri *>& nexttris() const noexcept
    { return pNextTris; }
    inline uint countNextTris() const noexcept
    { return pNextTris.size(); }


    inline std::vector<smtos::KProc *>::const_iterator kprocBegin() const noexcept
    { return pKProcs.begin(); }
    inline std::vector<smtos::KProc *>::const_iterator kprocEnd() const noexcept
    { return pKProcs.end(); }
    inline uint getStartKProcIdx() const noexcept
    {return startKProcIdx;}
    inline uint countKProcs() const noexcept
    { return nKProcs; }
    inline std::vector<smtos::KProc *> & kprocs() noexcept
    { return pKProcs; }
    
    inline KProc * getKProc(uint lidx)
    {
        if (hostRank != myRank) return nullptr;
        AssertLog(lidx < pKProcs.size());
        return pKProcs[lidx];
    }

    smtos::Reac * reac(uint lidx) const;
    
    virtual double getPoolOccupancy(uint lidx) const;
    virtual double getLastUpdate(uint lidx) const;
    virtual void resetPoolOccupancy();

    ////////////////////////////////////////////////////////////////////////
    
    // MPISTEPS
    bool getInHost() const;
    inline int getHost() const { return hostRank; }
    void setHost(int host, int rank);
    //void addSyncHost(int host);
    void setSolver(steps::mpi::tetopsplit::TetOpSplitP* solver);
    steps::mpi::tetopsplit::TetOpSplitP* solver() const;
    
    //virtual void sendSyncPools();
    //void recvSyncPools(int source);
    // MPISTEPS
    
    /////////////////////////// Dependency ////////////////////////////////
    // setup dependence for KProcs in this subvolume
    virtual void setupDeps();
    
    // check if kp_lidx in this vol depends on spec_gidx in WMVol kp_container
    virtual bool KProcDepSpecTet(uint kp_lidx, WmVol* kp_container, uint spec_gidx);
    // check if kp_lidx in this vol depends on spec_gidx in Tri kp_container
    virtual bool KProcDepSpecTri(uint kp_lidx, Tri* kp_container, uint spec_gidx);
    
    virtual void repartition(smtos::TetOpSplitP * tex, int rank, int host_rank);
    
protected:

    /// Use to store inprocess KProcs.
    std::vector<smtos::KProc *>          pKProcs;
    // Use to store indices of outporcess KProcs
    uint                                 startKProcIdx;
    uint                                 nKProcs;
    // The connected patch triangles.
    // Could be any number from zero to no upper limit- if this object is used
    // to descirbe a well-mixed compartment this may be a big number
    std::vector<smtos::Tri * >				pNextTris;


    ////////////////////////////////////////////////////////////////////////

    tetrahedron_id_t                    pIdx;

    steps::solver::Compdef            * pCompdef;

    double                              pVol;

    /// Numbers of molecules -- stored as uint.
    uint                              * pPoolCount{nullptr};
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags{nullptr};

    ///////// MPI STUFFS ////////////////////////////////////////////////////
    int                                 myRank;
    int                                 hostRank;
    steps::mpi::tetopsplit::TetOpSplitP         * pSol;
    
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_WMVOL_HPP

// END



