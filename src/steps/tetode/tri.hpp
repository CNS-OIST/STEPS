/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
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

// logging
#include "util/common.h"
#include <easylogging++.h>
#include "util/error.hpp"

// STEPS headers.
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
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

    Tri(triangle_id_t idx, steps::solver::Patchdef *patchdef, double area,
        double l0, double l1, double l2, double d0, double d1, double d2,
        tetrahedron_id_t tetinner, tetrahedron_id_t tetouter,
        triangle_id_t tri0, triangle_id_t tri1, triangle_id_t tri2);
    ~Tri();

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

    inline steps::solver::Patchdef * patchdef() const noexcept
    { return pPatchdef; }

    inline triangle_id_t idx() const noexcept
    { return pIdx; }

    inline double area() const noexcept
    { return pArea; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline stode::Tet * iTet() const noexcept
    { return pInnerTet; }

    inline stode::Tet * oTet() const noexcept
    { return pOuterTet; }

    inline stode::Tri * nextTri(uint i) const
    {
        AssertLog(i < 3);
        return pNextTri[i];
    }

    inline triangle_id_t tri(uint t) noexcept
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

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    double getOhmicI(double v,  stode::TetODE * solver) const;

    double getGHKI(double v,double dt, steps::tetode::TetODE * solver) const;

    /*
    inline uint * pools() const
    { return pPoolCount; }
    void setCount(uint lidx, uint count);


    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    */

    ////////////////////////////////////////////////////////////////////////

    /*
    inline std::vector<stex::KProc *>::const_iterator kprocBegin() const
    { return pKProcs.begin(); }
    inline std::vector<stex::KProc *>::const_iterator kprocEnd() const
    { return pKProcs.end(); }
    inline uint countKProcs() const
    { return pKProcs.size(); }

    stex::SReac * sreac(uint lidx) const;
    stex::VDepTrans * vdeptrans(uint lidx) const;
    stex::VDepSReac * vdepsreac(uint lidx) const;
    stex::GHKcurr * ghkcurr(uint lidx) const;

     */

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    triangle_id_t                       pIdx;

    steps::solver::Patchdef           * pPatchdef;

    /// Pointers to neighbouring tetrahedra.
    stode::Tet                       * pInnerTet{nullptr};
    stode::Tet                       * pOuterTet{nullptr};

    // Indices of two neighbouring tets; -1 if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    tetrahedron_id_t                   pTets[2];

    // Indices of neighbouring triangles.
    triangle_id_t                      pTris[3];

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
