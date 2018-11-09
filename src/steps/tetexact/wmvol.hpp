/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_TETEXACT_WMVOL_HPP
#define STEPS_TETEXACT_WMVOL_HPP 1

// STL headers.
#include <cassert>
#include <vector>
#include <fstream>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/compdef.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/solver/types.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace tetexact{

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class WmVol;
class Tri;
class Reac;
class Tetexact;

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
        uint idx, steps::solver::Compdef * cdef, double vol
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
    virtual void setupKProcs(stex::Tetexact * tex);

    virtual void setNextTri(stex::Tri *t);

    ////////////////////////////////////////////////////////////////////////

    virtual void reset();

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    inline steps::solver::Compdef * compdef() const
    { return pCompdef; }

    inline uint idx() const
    { return pIdx; }

    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume.
    ///
    inline double vol() const
    { return pVol; }


    ////////////////////////////////////////////////////////////////////////

    inline uint * pools() const
    { return pPoolCount; }
    void setCount(uint lidx, uint count);
    void incCount(uint lidx, int inc);

    // The concentration of species global index gidx in MOL PER l
    double conc(uint gidx) const;

    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<stex::Tri *>::const_iterator nexttriBegin() const
    { return pNextTris.begin(); }
    inline std::vector<stex::Tri *>::const_iterator nexttriEnd() const
    { return pNextTris.end(); }
    inline uint countNextTris() const
    { return pNextTris.size(); }
    inline const std::vector<stex::Tri *> & nexttris() const
    { return pNextTris; }


    inline std::vector<stex::KProc *>::const_iterator kprocBegin() const
    { return pKProcs.begin(); }
    inline std::vector<stex::KProc *>::const_iterator kprocEnd() const
    { return pKProcs.end(); }
    inline uint countKProcs() const
    { return pKProcs.size(); }
    inline std::vector<stex::KProc *> & kprocs()
    { return pKProcs; }

    stex::Reac * reac(uint lidx) const;

    ////////////////////////////////////////////////////////////////////////

protected:

    /// The kinetic processes.
    std::vector<stex::KProc *>          pKProcs;

    // The connected patch triangles.
    // Could be any number from zero to no upper limit- if this object is used
    // to descirbe a well-mixed compartment this may be a big number
    std::vector<stex::Tri * >                pNextTris;

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                 pIdx;

    steps::solver::Compdef            * pCompdef;

    double                              pVol;

    /// Numbers of molecules -- stored as uint.
    uint                              * pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags;

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_WMVOL_HPP

// END



