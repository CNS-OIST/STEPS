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

#include <iosfwd>
#include <string>

#include "fwd.hpp"
#include "model/fwd.hpp"

namespace steps::solver {

class GHKcurrdef {
  public:
    /// Constructor
    ///
    /// \param sd Defined state of the solver.
    /// \param gidx Global index of the GHK current.
    /// \param ghk Reference to the GHKcurr object.
    GHKcurrdef(Statedef& sd, ghkcurr_global_id gidx, model::GHKcurr& ghk);

    GHKcurrdef(const GHKcurrdef&) = delete;
    GHKcurrdef& operator=(const GHKcurrdef&) = delete;

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////

    /// Setup the object.
    void setup(const Statedef& sd);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: GHK CURRENT
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this ghk current.
    inline ghkcurr_global_id gidx() const noexcept {
        return pIdx;
    }

    /// Return the name of the ghk current.
    inline std::string const name() const noexcept {
        return pName;
    }

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: ION AND CHANNEL STATE
    ////////////////////////////////////////////////////////////////////////

    // Return the global index of the ion.
    spec_global_id ion() const;

    // Return real flux flag
    inline bool realflux() const noexcept {
        return pRealFlux;
    }

    // Return virtual outer concentration
    inline double voconc() const noexcept {
        return pVirtual_oconc;
    }

    // Return voltage-shift
    inline double vshift() const noexcept {
        return pVshift;
    }

    /// Return the calculated single channel permeability.
    inline double perm() const noexcept {
        return pPerm;
    }

    // Return the valence of the ion
    inline int valence() const noexcept {
        return pValence;
    }

    // Return the global index of the channel state
    spec_global_id chanstate() const;

    // For channel state
    int dep(spec_global_id gidx) const;
    bool req(spec_global_id gidx) const;

    // For species in inner and outer volumes
    int dep_v(spec_global_id gidx) const;
    bool req_v(spec_global_id gidx) const;

  private:
    /// The global index of this ghk current.
    const ghkcurr_global_id pIdx;

    /// The string identifier of this ghk current.
    const std::string pName;

    /// True if setup() has been called.
    bool pSetupdone{false};

    ////////////////////////////////////////////////////////////////////////
    // DATA: PARAMETERS
    ////////////////////////////////////////////////////////////////////////

    /// The channel state stored as a string, rather than ChanState pointer
    const std::string pChanState;

    /// The Ion stored as a string, rather than Spec pointer
    const std::string pIon;

    /// Flag whether the current is modelled as real movement of ions or not
    const bool pRealFlux;

    /// The virtual outer conc: if this is positive then the
    /// concentration in the outer compartment (if it exists) will be ignored
    const double pVirtual_oconc;

    /// The voltage-shift for the current calculation, defaults to zero.
    const double pVshift;

    /// The single-channel permeability  information.
    /// This is calculated internally from the conductance information
    /// supplied by user to GHKcurr object
    const double pPerm;

    /// The ion valence, copied for calculation of the GHK flux
    const int pValence;

    util::strongid_vector<spec_global_id, int> pSpec_DEP;

    util::strongid_vector<spec_global_id, int> pSpec_VOL_DEP;

    /// Global index of the channel state
    spec_global_id pSpec_CHANSTATE{};

    /// Global index of the ion.
    spec_global_id pSpec_ION{};
};

}  // namespace steps::solver
