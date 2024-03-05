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

#include "model/complex.hpp"
#include "statedef.hpp"
#include "util/common.hpp"

namespace steps::solver {

/// Defined Complex
class Complexdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the complex.
    /// \param d Pointer to the associated Complex object.
    Complexdef(complex_global_id idx, steps::model::Complex& d);

    /// Destructor
    ~Complexdef();

    ////////////////////////////////////////////////////////////////////////
    // CHECKPOINTING
    ////////////////////////////////////////////////////////////////////////
    /// checkpoint data
    void checkpoint(std::fstream& cp_file) const;

    /// restore data
    void restore(std::fstream& cp_file);

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this complex.
    complex_global_id gidx() const noexcept {
        return pIdx;
    }

    uint nbSubStates() const noexcept {
        return pNbSubStates;
    }

    /// Return the name of the complex.
    const std::string& name() const {
        return pName;
    }

    ////////////////////////////////////////////////////////////////////////
    // SOLVER METHODS: SETUP
    ////////////////////////////////////////////////////////////////////////
    /// Setup the object.
    ///
    /// This method is included for consistency with other def objects,
    /// but currently does nothing.
    void setup();

  private:
    const complex_global_id pIdx;
    const uint pNbSubStates;
    const std::string pName;
    bool pSetupdone{false};
};

}  // namespace steps::solver
