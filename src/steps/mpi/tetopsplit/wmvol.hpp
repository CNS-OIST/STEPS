/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2017 Okinawa Institute of Science and Technology, Japan.
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
        uint idx, steps::solver::Compdef * cdef, double vol, int rank, int host_rank
    );

    virtual ~WmVol(void);

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

    virtual void reset(void);

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    inline steps::solver::Compdef * compdef(void) const
    { return pCompdef; }

    inline uint idx(void) const
    { return pIdx; }

    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume.
    ///
    inline double vol(void) const
    { return pVol; }


    ////////////////////////////////////////////////////////////////////////

    inline uint * pools(void) const
    { return pPoolCount; }
    virtual void setCount(uint lidx, uint count, double period = 0.0);
    virtual void incCount(uint lidx, int inc, double period = 0.0);

    // The concentration of species global index gidx in MOL PER l
    double conc(uint gidx) const;

    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<smtos::Tri *>::const_iterator nexttriBegin(void) const
    { return pNextTris.begin(); }
    inline std::vector<smtos::Tri *>::const_iterator nexttriEnd(void) const
    { return pNextTris.end(); }
    inline uint countNextTris(void) const
    { return pNextTris.size(); }


    inline std::vector<smtos::KProc *>::const_iterator kprocBegin(void) const
    { return pKProcs.begin(); }
    inline std::vector<smtos::KProc *>::const_iterator kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint getStartKProcIdx(void)
    {return startKProcIdx;}
    inline uint countKProcs(void) const
    {
        return nKProcs;
    }
    inline std::vector<smtos::KProc *> & kprocs(void)
    { return pKProcs; }
    
    inline KProc * getKProc(uint lidx)
    {
        if (hostRank != myRank) return NULL;
        assert(lidx < pKProcs.size());
        return pKProcs[lidx];
    }

    smtos::Reac * reac(uint lidx) const;
    
    virtual double getPoolOccupancy(uint lidx);
    virtual double getLastUpdate(uint lidx);
    virtual void resetPoolOccupancy(void);

    ////////////////////////////////////////////////////////////////////////
    
    // MPISTEPS
    bool getInHost(void);
    int getHost(void) {return hostRank;}
    void setHost(int host, int rank);
    //void addSyncHost(int host);
    void setSolver(steps::mpi::tetopsplit::TetOpSplitP* solver);
    steps::mpi::tetopsplit::TetOpSplitP* solver(void);
    
    //virtual void sendSyncPools(void);
    //void recvSyncPools(int source);
    // MPISTEPS
    
    /////////////////////////// Dependency ////////////////////////////////
    // setup dependence for KProcs in this subvolume
    void setupDeps(void);
    
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

    uint                                 pIdx;

    steps::solver::Compdef            * pCompdef;

    double                              pVol;

    /// Numbers of molecules -- stored as uint.
    uint                              * pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags;

    ///////// MPI STUFFS ////////////////////////////////////////////////////
    int                                 hostRank;
    int                                 myRank;
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



