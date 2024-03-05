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
#include <fstream>
#include <vector>

// STEPS headers.
#include "kproc.hpp"
// #include "solver/compdef.hpp"
#include "util/vocabulary.hpp"

namespace steps::tetexact {

////////////////////////////////////////////////////////////////////////////////

// Forward declarations
class WmVol;
class Tri;
class Reac;
class Tetexact;

// Auxiliary declarations.
typedef WmVol* WmVolP;
typedef std::vector<WmVolP> WmVolPVec;
typedef WmVolPVec::iterator WmVolPVecI;
typedef WmVolPVec::const_iterator WmVolPVecCI;

////////////////////////////////////////////////////////////////////////////////

// Base class for the tetrahedrons in the mesh. This allows for compartments to
// be described as a well-mixed volume or comprised of tetrahedrons in the
// reaction-diffusion solver. Of course, if a compartment is well-mixed,
// any diffusion rules are ignored.
//
class WmVol {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    WmVol(tetrahedron_global_id idx, solver::Compdef* cdef, double vol);

    virtual ~WmVol();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    virtual void checkpoint(std::fstream& cp_file);

    /// restore data
    virtual void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Create the kinetic processes -- to be called when all tetrahedrons
    /// and triangles have been fully declared and connected.
    ///
    virtual void setupKProcs(Tetexact* tex);

    virtual void setNextTri(Tri* t);

    ////////////////////////////////////////////////////////////////////////

    virtual void reset();

    ////////////////////////////////////////////////////////////////////////
    // GENERAL INFORMATION
    ////////////////////////////////////////////////////////////////////////
    inline solver::Compdef* compdef() const noexcept {
        return pCompdef;
    }

    inline tetrahedron_global_id idx() const noexcept {
        return pIdx;
    }


    ////////////////////////////////////////////////////////////////////////
    // SHAPE & CONNECTIVITY INFORMATION.
    ////////////////////////////////////////////////////////////////////////

    /// Get the volume.
    ///
    inline double vol() const noexcept {
        return pVol;
    }


    ////////////////////////////////////////////////////////////////////////

    inline const auto& pools() const noexcept {
        return pPoolCount;
    }
    void setCount(solver::spec_local_id lidx, uint count);
    void incCount(solver::spec_local_id lidx, int inc);

    // The concentration of species global index gidx in MOL PER l
    double conc(solver::spec_global_id gidx) const;

    static const uint CLAMPED = 1;

    inline bool clamped(solver::spec_local_id lidx) const noexcept {
        return pPoolFlags[lidx] & CLAMPED;
    }
    void setClamped(solver::spec_local_id lidx, bool clamp);

    ////////////////////////////////////////////////////////////////////////

    inline uint countNextTris() const noexcept {
        return static_cast<uint>(pNextTris.size());
    }
    inline const std::vector<Tri*>& nexttris() const noexcept {
        return pNextTris;
    }

    inline uint countKProcs() const noexcept {
        return static_cast<uint>(pKProcs.size());
    }
    inline std::vector<KProc*>& kprocs() noexcept {
        return pKProcs;
    }

    Reac& reac(solver::reac_local_id lidx) const;

    ////////////////////////////////////////////////////////////////////////

  protected:
    /// The kinetic processes.
    std::vector<KProc*> pKProcs;

    // The connected patch triangles.
    // Could be any number from zero to no upper limit- if this object is used
    // to descirbe a well-mixed compartment this may be a big number
    std::vector<Tri*> pNextTris;

  private:
    ////////////////////////////////////////////////////////////////////////

    tetrahedron_global_id pIdx;

    solver::Compdef* pCompdef;

    double pVol;

    /// Numbers of molecules -- stored as uint.
    util::strongid_vector<solver::spec_local_id, uint> pPoolCount;
    /// Flags on these pools -- stored as machine word flags.
    util::strongid_vector<solver::spec_local_id, uint> pPoolFlags;
};

}  // namespace steps::tetexact
