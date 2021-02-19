/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2021 Okinawa Institute of Science and Technology, Japan.
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
#include <iostream>
#include <string>

// STEPS headers.
#include "steps/common.h"
#include "steps/error.hpp"
#include "steps/rng/create.hpp"
#include "steps/rng/mt19937.hpp"
#include "steps/rng/r123.hpp"

// logging
#include "easylogging++.h"
////////////////////////////////////////////////////////////////////////////////

namespace steps {
namespace rng {

RNGptr create(std::string rng_name, uint bufsize) {
    if (rng_name == "mt19937") {
        return RNGptr(new MT19937(bufsize));
    }
    else if (rng_name == "r123") {
        return RNGptr(new R123(bufsize));
    }
    else {
        ArgErrLog("Random number generator " + rng_name + " currently not included in STEPS.");
    }
}

RNGptr create_mt19937(uint bufsize) {
    return RNGptr(new MT19937(bufsize));
}

}  // namespace rng
}  // namespace steps

////////////////////////////////////////////////////////////////////////////////

// END
