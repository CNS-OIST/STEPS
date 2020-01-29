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


#ifndef STEPS_TETEXACT_TRI_HPP
#define STEPS_TETEXACT_TRI_HPP 1

// STL headers.
#include <cassert>
#include <vector>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "steps/common.h"
#include "steps/solver/patchdef.hpp"
#include "steps/tetexact/kproc.hpp"
#include "steps/solver/types.hpp"
////////////////////////////////////////////////////////////////////////////////

namespace steps{
namespace tetexact{

////////////////////////////////////////////////////////////////////////////////

namespace stex = steps::tetexact;

////////////////////////////////////////////////////////////////////////////////

// Forward declarations

class Tet;
class WmVol;
class Tri;
class SReac;
class SDiff;
class Tetexact;
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

    Tri(triangle_id_t idx, steps::solver::Patchdef * patchdef, double area,
        double l0, double l1, double l2, double d0, double d1, double d2,
        tetrahedron_id_t tetinner, tetrahedron_id_t tetouter,
        triangle_id_t tri0, triangle_id_t tri1, triangle_id_t tri2);
    ~Tri();

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
    void setInnerTet(stex::WmVol * t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(stex::WmVol * t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, stex::Tri * t);


    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(stex::Tetexact * tex, bool efield = false);

    /// Set all pool flags and molecular populations to zero.
    void reset();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GENERAL
    ////////////////////////////////////////////////////////////////////////

    inline steps::solver::Patchdef * patchdef() const
    { return pPatchdef; }

    inline triangle_id_t idx() const
    { return pIdx; }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline double area() const
    { return pArea; }

    inline stex::WmVol * iTet() const
    { return pInnerTet; }

    inline stex::WmVol * oTet() const
    { return pOuterTet; }

    inline stex::Tri * nextTri(uint i) const
    {
        AssertLog(i < 3);
        return pNextTri[i];
    }

    inline triangle_id_t tri(uint t)
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

    inline tetrahedron_id_t tet(uint t) const
    { return pTets[t]; }
    
    /// Find the direction index towards a neighbor triangle.
    ///
    int getTriDirection(triangle_id_t tidx);

    ////////////////////////////////////////////////////////////////////////

    // Set whether a direction is a diffusion boundary
    void setSDiffBndDirection(uint i);

    inline bool getSDiffBndDirection(uint idx) const
    { return pSDiffBndDirection[idx]; }

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

    inline uint * pools() const noexcept
    { return pPoolCount; }
    void setCount(uint lidx, uint count);
    void incCount(uint lidx, int inc);


    static const uint CLAMPED = 1;

    inline bool clamped(uint lidx) const noexcept
    { return pPoolFlags[lidx] & CLAMPED; }
    void setClamped(uint lidx, bool clamp);

    // Set a channel state relating to an ohmic current change.
    // 0th argument is oc local index, 1st argument is the local index
    // of the related channel state
    void setOCchange(uint oclidx, uint slidx, double dt, double simtime);

    ////////////////////////////////////////////////////////////////////////

    inline std::vector<stex::KProc *>::const_iterator kprocBegin() const noexcept
    { return pKProcs.begin(); }
    inline std::vector<stex::KProc *>::const_iterator kprocEnd() const noexcept
    { return pKProcs.end(); }
    inline std::vector<stex::KProc *> & kprocs() noexcept
    { return pKProcs; }
    inline const std::vector<stex::KProc *> & kprocs() const noexcept
    { return pKProcs; }
    inline uint countKProcs() const noexcept
    { return pKProcs.size(); }

    stex::SReac * sreac(uint lidx) const;
    stex::SDiff * sdiff(uint lidx) const;
    stex::VDepTrans * vdeptrans(uint lidx) const;
    stex::VDepSReac * vdepsreac(uint lidx) const;
    stex::GHKcurr * ghkcurr(uint lidx) const;

    ////////////////////////////////////////////////////////////////////////

private:

    ////////////////////////////////////////////////////////////////////////

    triangle_id_t                                 pIdx;

    steps::solver::Patchdef           * pPatchdef;

    /// Pointers to neighbouring tetrahedra.
    stex::WmVol                       * pInnerTet;
    stex::WmVol                       * pOuterTet;

    // Indices of two neighbouring tets; -1 if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    tetrahedron_id_t                    pTets[2];

    // Indices of neighbouring triangles.
    triangle_id_t                       pTris[3];

    /// POinters to neighbouring triangles
    stex::Tri                         * pNextTri[3];


    double                              pArea;

    // Neighbour information- needed for surface diffusion
    double                              pLengths[3];
    double                              pDist[3];

    bool                                pSDiffBndDirection[3];

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

    ////////////////////////////////////////////////////////////////////////

};

////////////////////////////////////////////////////////////////////////////////

}
}

////////////////////////////////////////////////////////////////////////////////

#endif

// STEPS_TETEXACT_TRI_HPP

// END
