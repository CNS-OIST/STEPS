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


#ifndef STEPS_TETEXACT_KPROC_HPP
#define STEPS_TETEXACT_KPROC_HPP 1

////////////////////////////////////////////////////////////////////////////////


// STL headers.
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/types.hpp"
#include "steps/rng/rng.hpp"
//#include "tetexact.hpp"

// Tetexact CR header
#include "steps/tetexact/crstruct.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace tetexact{

////////////////////////////////////////////////////////////////////////////////

//Forward declaration
class Tet;
class Tri;
class WmVol;
class KProc;
class Tetexact;

////////////////////////////////////////////////////////////////////////////////

typedef KProc *                         KProcP;
typedef std::vector<KProcP>             KProcPVec;
typedef KProcPVec::iterator             KProcPVecI;
typedef KProcPVec::const_iterator       KProcPVecCI;

////////////////////////////////////////////////////////////////////////////////

class KProc

{
public:

    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    KProc(void);
    virtual ~KProc(void);

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

    inline bool active(void) const
    { return !(pFlags & INACTIVATED); }
    inline bool inactive(void) const
    { return (pFlags & INACTIVATED); }
    void setActive(bool active);

    inline uint flags(void) const
    { return pFlags; }

    ////////////////////////////////////////////////////////////////////////

    uint schedIDX(void) const
    { return pSchedIDX; }

    void setSchedIDX(uint idx)
    { pSchedIDX = idx; }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    /// This function is called when all kproc objects have been created,
    /// allowing the kproc to pre-compute its SchedIDXVec.
    ///
    virtual void setupDeps(void) = 0;

    virtual bool depSpecTet(uint gidx, steps::tetexact::WmVol * tet) = 0;
    virtual bool depSpecTri(uint gidx, steps::tetexact::Tri * tri) = 0;

    /// Reset this Kproc.
    ///
    virtual void reset(void) = 0;

    // Recompute the Ccst for this KProc
    virtual void resetCcst(void) const;

    /// Compute the rate for this kproc (its propensity value).
    ///
    virtual double rate(steps::tetexact::Tetexact * solver = 0)  = 0;

    // Return the ccst for this kproc
    // NOTE: not pure for this solver because doesn't make sense for Diff
    virtual double c(void) const;

    // Return the h value for this kproc (number of available reaction channels)
    // NOTE: not pure for this solver because doesn;t make sense for Diff
    virtual double h(void);

    /// Apply a single discrete instance of the kinetic process, returning
    /// a vector of kproc schedule indices that need to be updated as a
    /// result.
    ///
    // NOTE: Random number generator available to this function for use
    // by Diff
    virtual std::vector<KProc*> const & apply(steps::rng::RNG * rng, double dt, double simtime) = 0;

    virtual uint updVecSize(void) const = 0;

    ////////////////////////////////////////////////////////////////////////

    uint getExtent(void) const;
    void resetExtent(void);

    ////////////////////////////////////////////////////////////////////////
    /*
    // Return a pointer to the corresponding Reacdef Diffdef or SReacdef
    // object
    // Separate methods to avoid making a base KProcdef class
    //
    virtual steps::solver::Reacdef * defr(void) const;
    virtual steps::solver::SReacdef * defsr(void) const;
    */// compileerror; // check this

    ////////////////////////////////////////////////////////////////////////

    // data for CR SSA
    CRKProcData                         crData;

protected:

    uint                                rExtent;

    ////////////////////////////////////////////////////////////////////////

    uint                                pFlags;

    uint                                pSchedIDX;

    ////////////////////////////////////////////////////////////////////////
};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_KPROC_HPP
