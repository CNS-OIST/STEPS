/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2026 Okinawa Institute of Science and Technology, Japan.
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
#include <string>

// STEPS headers.
#include "fwd.hpp"
#include "model/complex.hpp"
#include "util/common.hpp"

////////////////////////////////////////////////////////////////////////////////

namespace steps::dist {

/// Defined Complex
class Complexdef {
  public:
    /// Constructor
    ///
    /// \param sd State of the solver.
    /// \param idx Global index of the complex.
    /// \param d Reference to the associated Complex object.
    Complexdef(const Statedef& sd, model::complex_id idx, const steps::model::Complex& d);

    /// Reset the def-object values to model defaults
    void reset() {}

    /// Destructor
    ~Complexdef();

    ////////////////////////////////////////////////////////////////////////
    // DATA ACCESS: COMPLEX
    ////////////////////////////////////////////////////////////////////////

    /// Return the global index of this complex.
    inline model::complex_id idx() const noexcept {
        return pIdx;
    }

    inline uint nbSubStates() const noexcept {
        return pNbSubStates;
    }

    ////////////////////////////////////////////////////////////////////////

  private:
    ////////////////////////////////////////////////////////////////////////

    model::complex_id pIdx;
    uint pNbSubStates;

    ////////////////////////////////////////////////////////////////////////
};

}  // namespace steps::dist
