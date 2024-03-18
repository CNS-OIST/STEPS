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
#include <fstream>
#include <set>
#include <vector>

// STEPS headers.
#include "crstruct.hpp"
#include "rng/rng.hpp"
#include "solver/fwd.hpp"
#include "util/collections.hpp"

namespace steps::tetexact {

// Forward declaration
class Tet;
class Tri;
class WmVol;
class KProc;
class Tetexact;

typedef KProc* KProcP;
typedef std::vector<KProcP> KProcPVec;

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
    virtual void checkpoint(std::fstream& cp_file);

    /// restore data
    virtual void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS
    ////////////////////////////////////////////////////////////////////////

    static const int INACTIVATED = 1;

    inline bool active() const noexcept {
        return !(pFlags & INACTIVATED);
    }
    inline bool inactive() const noexcept {
        return static_cast<bool>(pFlags & INACTIVATED);
    }
    void setActive(bool active);

    inline uint flags() const {
        return pFlags;
    }

    ////////////////////////////////////////////////////////////////////////

    solver::kproc_global_id schedIDX() const {
        return pSchedIDX;
    }

    void setSchedIDX(solver::kproc_global_id idx) {
        pSchedIDX = idx;
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    /// This function is called when all kproc objects have been created,
    /// allowing the kproc to pre-compute its SchedIDXVec.
    ///
    virtual void setupDeps() = 0;

    virtual bool depSpecTet(solver::spec_global_id gidx, steps::tetexact::WmVol* tet) = 0;
    virtual bool depSpecTri(solver::spec_global_id gidx, steps::tetexact::Tri* tri) = 0;

    /// Reset this Kproc.
    ///
    virtual void reset() = 0;

    // Recompute the Ccst for this KProc
    virtual void resetCcst() const;

    /// Compute the rate for this kproc (its propensity value).
    ///
    virtual double rate(steps::tetexact::Tetexact* solver = nullptr) = 0;

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
    virtual std::vector<KProc*> const& apply(const rng::RNGptr& rng, double dt, double simtime) = 0;

    virtual uint updVecSize() const = 0;

    ////////////////////////////////////////////////////////////////////////

    unsigned long long getExtent() const;
    void resetExtent();

    ////////////////////////////////////////////////////////////////////////

    // data for CR SSA
    CRKProcData crData;

  protected:
    unsigned long long rExtent{0};

    ////////////////////////////////////////////////////////////////////////

    uint pFlags{0};

    solver::kproc_global_id pSchedIDX{};

    ////////////////////////////////////////////////////////////////////////
};

inline bool operator<(const KProc& lhs, const KProc& rhs) {
    return lhs.schedIDX() < rhs.schedIDX();
}

using KProcPSet = std::set<KProc*, util::DerefPtrLess<KProc>>;

}  // namespace steps::tetexact

namespace std {
// Compilation trap in case std::set<KProc*> is used in the code that will sort KProc instances
// by their pointer addresses, not their schedule identifier.
// Prefer KProcPSet for such usage
template <>
class set<steps::tetexact::KProc*> {};
}  // namespace std
