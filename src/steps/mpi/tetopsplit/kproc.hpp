/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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


/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#ifndef STEPS_MPI_TETOPSPLIT_KPROC_HPP
#define STEPS_MPI_TETOPSPLIT_KPROC_HPP 1

////////////////////////////////////////////////////////////////////////////////


// STL headers.
#include <vector>
#include <fstream>
#include <set>
#include <random>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/rng/rng.hpp"
//#include "tetopsplit.hpp"

// TetOpSplitP CR header
#include "steps/mpi/tetopsplit/crstruct.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace mpi{
namespace tetopsplit{

////////////////////////////////////////////////////////////////////////////////

//Forward declaration
class Tet;
class Tri;
class WmVol;
class KProc;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

typedef KProc *                         KProcP;
typedef std::vector<KProcP>             KProcPVec;
typedef KProcPVec::iterator             KProcPVecI;
typedef KProcPVec::const_iterator       KProcPVecCI;

typedef std::set<KProcP>                KProcPSet;
typedef KProcPSet::iterator             KProcPSetI;
typedef KProcPSet::const_iterator       KProcPSetCI;

////////////////////////////////////////////////////////////////////////////////

enum TYPE {KP_REAC, KP_SREAC, KP_DIFF, KP_SDIFF, KP_GHK, KP_VDEPSREAC, KP_VDEPTRANS};

class KProc

{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    KProc();
    virtual ~KProc();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    virtual void checkpoint(std::fstream & cp_file) = 0;

    /// restore data
    virtual void restore(std::fstream & cp_file) = 0;

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////


    static const int INACTIVATED = 1;

    inline bool active() const noexcept
    { return !(pFlags & INACTIVATED); }
    inline bool inactive() const noexcept
    { return (pFlags & INACTIVATED); }
    void setActive(bool active);

    inline uint flags() const noexcept
    { return pFlags; }

    ////////////////////////////////////////////////////////////////////////

    uint schedIDX() const noexcept
    { return pSchedIDX; }

    void setSchedIDX(uint idx) noexcept
    { pSchedIDX = idx; }
    
    uint getType() const noexcept { return type; }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    /// This function is called when all kproc objects have been created,
    /// allowing the kproc to pre-compute its SchedIDXVec.
    ///
    virtual void setupDeps() = 0;

    virtual bool depSpecTet(uint gidx, steps::mpi::tetopsplit::WmVol * tet = nullptr) = 0;
    virtual bool depSpecTri(uint gidx, steps::mpi::tetopsplit::Tri * tri = nullptr) = 0;

    /// Reset this Kproc.
    ///
    virtual void reset() = 0;

    // Recompute the Ccst for this KProc
    virtual void resetCcst();

    /// Compute the rate for this kproc (its propensity value).
    ///
    virtual double rate(steps::mpi::tetopsplit::TetOpSplitP * solver = nullptr)  = 0;
    virtual double getScaledDcst(steps::mpi::tetopsplit::TetOpSplitP * solver = nullptr) const = 0;

    // Return the ccst for this kproc
    // NOTE: not pure for this solver because doesn't make sense for Diff
    virtual double c() const;

    // Return the h value for this kproc (number of available reaction channels)
    // NOTE: not pure for this solver because doesn;t make sense for Diff
    virtual double h();

    /// Apply a single discrete instance of the kinetic process, returning
    /// a vector of kproc schedule indices that need to be updated as a
    /// result.
    ///
    // NOTE: Random number generator available to this function for use
    // by Diff

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
    virtual int apply(const rng::RNGptr &rng);
    virtual int apply(const rng::RNGptr &rng, uint nmolcs);
    virtual void apply(const rng::RNGptr &rng, double dt, double simtime, double period);
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
    virtual int apply(const rng::RNGptr &rng);
    virtual int apply(const rng::RNGptr &rng, uint nmolcs);
    virtual void apply(const rng::RNGptr &rng, double dt, double simtime, double period);
#pragma GCC diagnostic pop
#endif

    virtual std::vector<KProc*> const & getLocalUpdVec(int direction = -1) const;
    virtual std::vector<uint> const & getRemoteUpdVec(int direction = -1) const;

    // Intended for reactions within the SSA
    // Special case for SReacs where dt and simtime are needed if Ohmic Currents are involved, i.e. a
    // Surface reaction can open or close an ohmic current channel

    virtual void resetOccupancies();
    
    virtual bool getInHost() const = 0;
    virtual int getHost() const = 0;
	
	//virtual std::vector<KProc*> const & getSharedUpd();
    ////////////////////////////////////////////////////////////////////////

    unsigned long long getExtent() const;
    void resetExtent();

    ////////////////////////////////////////////////////////////////////////
    /*
    // Return a pointer to the corresponding Reacdef Diffdef or SReacdef
    // object
    // Separate methods to avoid making a base KProcdef class
    //
    virtual steps::solver::Reacdef * defr() const;
    virtual steps::solver::SReacdef * defsr() const;
    */// compileerror; // check this

    ////////////////////////////////////////////////////////////////////////

    // data for CR SSA
    CRKProcData                         crData;

protected:

    unsigned long long                  rExtent;

    ////////////////////////////////////////////////////////////////////////

    uint                                pFlags;

    uint                                pSchedIDX;
    
    uint                                 type;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

}
}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_MPI_TETOPSPLIT_KPROC_HPP
