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

#include <cassert>
#include <cmath>

// STEPS headers.
#include "math/point.hpp"
#include "math/sphere.hpp"

namespace steps::math {

position_rel_to_ves sphere_unit_randsurfpos(rng::RNGptr rng) {
    // Positions are relative to the centre of the sphere
    // Initial position is randomised
    double r1 = rng->getUnfII() * 2 - 1;
    double r2 = rng->getUnfII() * 2 - 1;

    // First part of the algorithm is to ensure the sum of the squares
    // are less than 1
    while ((r1 * r1) + (r2 * r2) >= 1.0) {
        r1 = rng->getUnfII() * 2 - 1;
        r2 = rng->getUnfII() * 2 - 1;
    }

    double x = 2 * r1 * std::sqrt(1 - (r1 * r1) - (r2 * r2));
    double y = 2 * r2 * std::sqrt(1 - (r1 * r1) - (r2 * r2));
    double z = 1 - 2 * ((r1 * r1) + (r2 * r2));

    position_rel_to_ves ran_loc{x, y, z};

    return ran_loc;
}

}  // namespace steps::math
