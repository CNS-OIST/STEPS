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

// STEPS headers.
#include "math/point.hpp"
#include "rng/rng.hpp"
#include "util/common.hpp"

namespace steps::math {

/** Generate random position on unit sphere tetrahedron barycentre.
 *
 * \param r1, r1 Uniform random numbers on -1, 1.
 * \return random point on surface.
 */
position_rel_to_ves sphere_unit_randsurfpos(rng::RNGptr rng);

}  // namespace steps::math