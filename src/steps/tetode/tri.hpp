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


#ifndef STEPS_TETODE_TRI_HPP
#define STEPS_TETODE_TRI_HPP 1

// STL headers.
#include <cassert>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/solver/patchdef.hpp"
#include "steps/solver/types.hpp"
//#include "../tetode/tetode.hpp"
// logging
#include "third_party/easyloggingpp/src/easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

 namespace steps {
 namespace tetode {

////////////////////////////////////////////////////////////////////////////////

namespace stode = steps::tetode;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations

class Tet;
class Tri;
class TetODE;

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

    Tri(uint idx, steps::solver::Patchdef * patchdef, double area,
            double l0, double l1, double l2, double d0, double d1, double d2,
            int tetinner, int tetouter, int tri0, int tri1, int tri2);
    ~Tri(void);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Set pointer to the 'inside' neighbouring tetrahedron.
    ///
    void setInnerTet(stode::Tet * t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(stode::Tet * t);

    /// Set pointer to the next neighbouring triangle.
     void setNextTri(uint i, stode::Tri * t);

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream & cp_file);

    /// restore data
    void restore(std::fstream & cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Patchdef * patchdef(void) const
    { return pPatchdef; }

    inline uint idx(void) const
    { return pIdx; }

    inline double area(void) const
    { return pArea; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline stode::Tet * iTet(void) const
    { return pInnerTet; }

    inline stode::Tet * oTet(void) const
    { return pOuterTet; }

    inline stode::Tri * nextTri(uint i) const
    {
        AssertLog(i < 3);
        return pNextTri[i];
    }

    inline int tri(uint t)
    { return pTris[t]; }

    /// Get the length of a boundary bar.
    ///
    inline double length(uint i) const
    { return pLengths[i]; }

    /// Get the distance to the centroid of the next neighbouring
    /// triangle.
    ///
    inline double dist(uint i) const
    { return pDist[i]; }

    inline int tet(uint t) const
    { return pTets[t]; }

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    double getOhmicI(double v,  stode::TetODE * solver) const;

    double getGHKI(double v,double dt, steps::tetode::TetODE * solver) const;

    /*
    inline uint * pools(void) const
    { return pPoolCount; }
    void setCount(uint lidx, uint count);


    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    */

    ////////////////////////////////////////////////////////////////////////

    /*
    inline std::vector<stex::KProc *>::const_iterator kprocBegin(void) const
    { return pKProcs.begin(); }
    inline std::vector<stex::KProc *>::const_iterator kprocEnd(void) const
    { return pKProcs.end(); }
    inline uint countKProcs(void) const
    { return pKProcs.size(); }

    stex::SReac * sreac(uint lidx) const;
    stex::VDepTrans * vdeptrans(uint lidx) const;
    stex::VDepSReac * vdepsreac(uint lidx) const;
    stex::GHKcurr * ghkcurr(uint lidx) const;

     */

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    uint                                 pIdx;

    steps::solver::Patchdef           * pPatchdef;

    /// Pointers to neighbouring tetrahedra.
    stode::Tet                       * pInnerTet;
    stode::Tet                       * pOuterTet;

    // Indices of two neighbouring tets; -1 if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    int                                 pTets[2];

    // Indices of neighbouring triangles.
    int                                 pTris[3];

    /// POinters to neighbouring triangles
    stode::Tri                        * pNextTri[3];

    double                              pArea;

    // Neighbour information- needed for surface diffusion
    double                              pLengths[3];
    double                              pDist[3];

    /*
    /// Numbers of molecules -- stored as machine word integers.
    uint                              * pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    uint                              * pPoolFlags;

    /// The kinetic processes.
    std::vector<stex::KProc *>          pKProcs;

    /// For the EFIELD calculation. An integer storing the amount of
    /// elementary charge from inner tet to outer tet (positive if
    /// net flux is positive, negative if net flux is negative) for
    /// one EField time-step.
    // NOTE: Now arrays so as to separate into different GHK currs,
    // for data access
    int                                  * pECharge;

    // to store the latest ECharge, so that the info is available to solver
    int                               * pECharge_last;


    // Store the Ohmic currents' channel's opening information by OC local indices
    // and the time since the related channel state changed it's number
    // The pOCchan_timeintg stores number of channel open * opening time
    // so at the end of the step this number/Efield dt will give the
    // mean number of channels open
    double                               * pOCchan_timeintg;

    double                               * pOCtime_upd;

     */

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_TRI_HPP

// END
