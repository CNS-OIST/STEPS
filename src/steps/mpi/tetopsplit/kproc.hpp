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

/*
 *  Last Changed Rev:  $Rev$
 *  Last Changed Date: $Date$
 *  Last Changed By:   $Author$
 */

#pragma once

// STL headers.
#include <fstream>
#include <functional>
#include <random>
#include <set>
#include <vector>

// STEPS headers.
#include "rng/rng.hpp"
#include "solver/fwd.hpp"
#include "util/collections.hpp"

// TetOpSplitP CR header
#include "crstruct.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::mpi::tetopsplit {

////////////////////////////////////////////////////////////////////////////////

// Forward declaration
class Tet;
class Tri;
class WmVol;
class KProc;
class TetOpSplitP;

////////////////////////////////////////////////////////////////////////////////

typedef KProc* KProcP;
typedef std::vector<KProcP> KProcPVec;

////////////////////////////////////////////////////////////////////////////////

enum TYPE {
    KP_REAC,
    KP_SREAC,
    KP_DIFF,
    KP_SDIFF,
    KP_GHK,
    KP_VDEPSREAC /*, KP_VDEPTRANS*/
};

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
        return (pFlags & INACTIVATED) == 0;
    }
    inline bool inactive() const noexcept {
        return (pFlags & INACTIVATED) != 0;
    }
    void setActive(bool active);

    inline uint flags() const noexcept {
        return pFlags;
    }

    ////////////////////////////////////////////////////////////////////////

    solver::kproc_global_id schedIDX() const noexcept {
        return pSchedIDX;
    }

    void setSchedIDX(solver::kproc_global_id idx) noexcept {
        pSchedIDX = idx;
    }

    uint getType() const noexcept {
        return type;
    }

    ////////////////////////////////////////////////////////////////////////
    // VIRTUAL INTERFACE METHODS
    ////////////////////////////////////////////////////////////////////////

    /// This function is called when all kproc objects have been created,
    /// allowing the kproc to pre-compute its SchedIDXVec.
    ///
    virtual void setupDeps() = 0;

    /// Reset this Kproc.
    ///
    virtual void reset() = 0;

    /*
    virtual bool depSpecTet(solver::spec_global_id gidx,
    mpi::tetopsplit::WmVol * tet = nullptr) = 0; virtual bool
    depSpecTri(solver::spec_global_id gidx, mpi::tetopsplit::Tri * tri =
    nullptr) = 0;
        */
    // Recompute the Ccst for this KProc
    virtual void resetCcst();

    /// Compute the rate for this kproc (its propensity value).
    ///
    virtual double rate(mpi::tetopsplit::TetOpSplitP* solver = nullptr) = 0;
    virtual double getScaledDcst(mpi::tetopsplit::TetOpSplitP* solver = nullptr) const = 0;

    // Return the ccst for this kproc
    // NOTE: not pure for this solver because doesn't make sense for Diff
    virtual double c() const;

    // Return the h value for this kproc (number of available reaction channels)
    // NOTE: not pure for this solver because doesn;t make sense for Diff
    virtual double h();

    unsigned long long getExtent() const;
    void resetExtent();

    ////////////////////////////////////////////////////////////////////////

    /// Apply a single discrete instance of the kinetic process, returning
    /// a vector of kproc schedule indices that need to be updated as a
    /// result.
    ///
    // NOTE: Random number generator available to this function for use
    // by Diff

#if defined(__clang__)
#pragma clang diagnostic push
#pragma clang diagnostic ignored "-Woverloaded-virtual"
    virtual int apply(const rng::RNGptr& rng);
    virtual int apply(const rng::RNGptr& rng, uint nmolcs);
    virtual void apply(const rng::RNGptr& rng, double dt, double simtime, double period);
#pragma clang diagnostic pop
#elif defined(__GNUC__) || defined(__GNUG__)
#pragma GCC diagnostic push
#pragma GCC diagnostic ignored "-Woverloaded-virtual"
    virtual int apply(const rng::RNGptr& rng);
    virtual int apply(const rng::RNGptr& rng, uint nmolcs);
    virtual void apply(const rng::RNGptr& rng, double dt, double simtime, double period);
#pragma GCC diagnostic pop
#endif

    ////////////////////// MPI FUNCTIONALITY /////////////////////

    virtual std::vector<KProc*> const& getLocalUpdVec(int direction = -1) const;
    virtual std::vector<solver::kproc_global_id> const& getRemoteUpdVec(int direction = -1) const;

    // Intended for reactions within the SSA
    // Special case for SReacs where dt and simtime are needed if Ohmic Currents
    // are involved, i.e. a Surface reaction can open or close an ohmic current
    // channel

    virtual void resetOccupancies();

    virtual bool getInHost() const = 0;
    virtual int getHost() const = 0;

    ////////////////////////////////////////////////////////////////////////

    // data for CR SSA
    CRKProcData crData;

  protected:
    unsigned long long rExtent;

    ////////////////////////////////////////////////////////////////////////

    uint pFlags;

    solver::kproc_global_id pSchedIDX;

    uint type{};

    ////////////////////////////////////////////////////////////////////////
};

inline bool operator<(const KProc& lhs, const KProc& rhs) {
    return lhs.schedIDX() < rhs.schedIDX();
}

using KProcPSet = std::set<KProc*, util::DerefPtrLess<KProc>>;

}  // namespace steps::mpi::tetopsplit

namespace std {
// Compilation trap in case std::set<KProc*> is used in the code that will sort KProc instances
// by their pointer addresses, not their schedule identifier.
// Prefer KProcPSet for such usage
template <>
class set<steps::mpi::tetopsplit::KProc*> {};
}  // namespace std
