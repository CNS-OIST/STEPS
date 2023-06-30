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

#include "rng/std_mt19937.hpp"

// logging
#include <easylogging++.h>

namespace steps::rng {

STDMT19937::STDMT19937(uint bufsize)
    : RNG(bufsize) {}

STDMT19937::~STDMT19937() = default;

void STDMT19937::concreteInitialize(ulong seed) {
    rng_.seed(seed);
}

void STDMT19937::concreteFillBuffer() {
    for (uint i = 0; i < rSize; ++i) {
        rBuffer[i] = rng_();
    }
}

////////////////////////////////////////////////////////////////////////////////

void STDMT19937::checkpoint(std::ostream& cp_file) const {
    RNG::checkpoint(cp_file);
    CLOG(INFO, "general_log") << "Warning - std_mt19937 checkpointing is not implemented, runs "
                                 "might not be reproducible "
                              << std::endl;
    // TODO Implement checkpoint for std_mt19937
}

////////////////////////////////////////////////////////////////////////////////

void STDMT19937::restore(std::istream& cp_file) {
    RNG::restore(cp_file);
    // TODO Implement checkpoint for std_mt19937
}

}  // namespace steps::rng
