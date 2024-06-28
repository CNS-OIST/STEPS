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
////////////////////////////////////////////////////////////////////////////////

// Standard library & STL headers.
#include <math.h>
#include <vector>

// STEPS headers.
#include "rng/rng.hpp"

namespace steps::mpi::tetvesicle {

class Qtable {
  public:
    ////////////////////////////////////////////////////////////////////////
    // OBJECT CONSTRUCTION & DESTRUCTION
    ////////////////////////////////////////////////////////////////////////

    Qtable(unsigned int size, double tau, const rng::RNGptr& r);

    /// checkpoint data
    void checkpoint(std::fstream& cp_file);

    /// restore data
    void restore(std::fstream& cp_file);

    void setup() noexcept;

    // Intended to be used if tablesize or tau changes during simulation
    void reinit(unsigned int size, double tau) noexcept;

    double getPhi() const noexcept;

    inline double getTau() const noexcept {
        return pTau;
    }


    ////////////////////////////////////////////////////////////////////////

  private:
    double Q(double theta) const noexcept;
    unsigned int pTablesize;
    double pTau;

    // The interpolation table
    std::vector<double> pX_interp;
    std::vector<double> pQ_values;

    // RNG stuff
    const rng::RNGptr rng;
};

}  // namespace steps::mpi::tetvesicle
