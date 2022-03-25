/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2022 Okinawa Institute of Science and Technology, Japan.
#    Copyright (C) 2003-2006 University of Antwerp, Belgium.
#    
#    See the file AUTHORS for details.
#    This file is part of STEPS.
#    
#    STEPS is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License version 2,
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


// Standard library & STL headers.
#include <memory>
#include <string>

// STEPS headers.
#include "util/common.h"
#include "rng.hpp"

namespace steps {
namespace rng {

/// Create a random number generator with name rng_name and return as RNG object.
///
/// \param rng_name Name of the random number generator.
/// \param buffsize Size of buffer.
RNGptr create(std::string rng_name, uint bufsize);

/// Create a MT19937 random number generator and return as RNG object.
///
/// \param buffsize Size of buffer.
RNGptr create_mt19937(uint bufsize);

} // namespace rng
} // namespace steps
