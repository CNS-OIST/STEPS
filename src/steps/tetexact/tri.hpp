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
#include <cassert>
#include <vector>

// logging
#include <easylogging++.h>

// STEPS headers.
#include "kproc.hpp"
#include "solver/patchdef.hpp"
#include "solver/types.hpp"
#include "util/common.hpp"

namespace steps::tetexact {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations

class Tet;
class WmVol;
class Tri;
class SReac;
class SDiff;
class Tetexact;
class VDepSReac;
class GHKcurr;

// Auxiliary declarations.
typedef Tri* TriP;
typedef std::vector<TriP> TriPVec;
typedef TriPVec::iterator TriPVecI;
typedef TriPVec::const_iterator TriPVecCI;

class Tri {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Tri(triangle_global_id idx,
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
        triangle_global_id tri2);
    ~Tri();

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
    void setInnerTet(WmVol* t);

    /// Set pointer to the 'outside' neighbouring tetrahedron.
    ///
    void setOuterTet(WmVol* t);

    /// Set pointer to the next neighbouring triangle.
    void setNextTri(uint i, Tri* t);


    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    void setupKProcs(Tetexact* tex, bool efield = false);

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

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: SHAPE & CONNECTIVITY
    ////////////////////////////////////////////////////////////////////////

    inline double area() const noexcept {
        return pArea;
    }

    inline WmVol* iTet() const noexcept {
        return pInnerTet;
    }

    inline WmVol* oTet() const noexcept {
        return pOuterTet;
    }

    inline Tri* nextTri(uint i) const {
        return pNextTri.at(i);
    }

    inline triangle_global_id tri(uint t) const {
        return pTris.at(t);
    }

    /// Get the length of a boundary bar.
    ///
    inline double length(uint i) const {
        return pLengths.at(i);
    }

    /// Get the distance to the centroid of the next neighbouring
    /// triangle.
    ///
    inline double dist(uint i) const {
        return pDist.at(i);
    }

    inline tetrahedron_global_id tet(uint t) const {
        return pTets.at(t);
    }

    /// Find the direction index towards a neighbor triangle.
    ///
    int getTriDirection(triangle_global_id tidx) const noexcept;

    ////////////////////////////////////////////////////////////////////////

    // Set whether a direction is a diffusion boundary
    void setSDiffBndDirection(uint i);

    inline bool getSDiffBndDirection(uint idx) const {
        return pSDiffBndDirection.at(idx);
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: EFIELD
    ////////////////////////////////////////////////////////////////////////

    // Local index of GHK current given
    void incECharge(solver::ghkcurr_local_id lidx, int charge);

    // Should be called at the beginning of every EField time-step
    void resetECharge(double dt, double efdt, double t);

    // reset the Ohmic current opening time integral info, also should be
    // called just before commencing or just after completing an EField dt
    void resetOCintegrals();

    double computeI(double v, double dt, double simtime, double efdt);

    double getOhmicI(double v, double dt) const;
    double getOhmicI(solver::ohmiccurr_local_id lidx, double v, double dt) const;

    double getGHKI() const;
    double getGHKI(solver::ghkcurr_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////
    // MAIN FUNCTIONALITY
    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }
    void setCount(solver::spec_local_id lidx, uint count);
    void incCount(solver::spec_local_id lidx, int inc);


    static const uint CLAMPED = 1;

    inline bool clamped(solver::spec_local_id lidx) const noexcept {
        return (pPoolFlags[lidx] & CLAMPED) != 0;
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

    inline std::vector<KProc*>& kprocs() noexcept {
        return pKProcs;
    }
    inline const std::vector<KProc*>& kprocs() const noexcept {
        return pKProcs;
    }
    inline uint countKProcs() const noexcept {
        return pKProcs.size();
    }

    SReac& sreac(solver::sreac_local_id lidx) const;
    SDiff& sdiff(solver::surfdiff_local_id lidx) const;
    VDepSReac& vdepsreac(solver::vdepsreac_local_id lidx) const;
    GHKcurr& ghkcurr(solver::ghkcurr_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    triangle_global_id pIdx;

    solver::Patchdef* pPatchdef;

    /// Pointers to neighbouring tetrahedra.
    WmVol* pInnerTet{nullptr};
    WmVol* pOuterTet{nullptr};

    // Indices of two neighbouring tets; -1 if surface triangle (if
    // triangle's patch is on the surface of the mesh, quite often the case)
    std::array<tetrahedron_global_id, 2> pTets;

    // Indices of neighbouring triangles.
    std::array<triangle_global_id, 3> pTris;

    /// Pointers to neighbouring triangles
    std::array<Tri*, 3> pNextTri;


    double pArea;

    // Neighbour information- needed for surface diffusion
    std::array<double, 3> pLengths;
    std::array<double, 3> pDist;

    std::array<bool, 3> pSDiffBndDirection;

    /// Numbers of molecules -- stored as machine word integers.
    util::strongid_vector<solver::spec_local_id, uint> pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    util::strongid_vector<solver::spec_local_id, uint> pPoolFlags;

    /// The kinetic processes.
    std::vector<KProc*> pKProcs;

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

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::tetexact
