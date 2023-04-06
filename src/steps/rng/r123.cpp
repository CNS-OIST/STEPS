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


// Standard library & STL headers.
#include <cassert>
#include <cstdint>
#include <sstream>
#include <string>

// STEPS headers.
#include "r123.hpp"
// util
#include "util/error.hpp"
// logging
#include <easylogging++.h>
////////////////////////////////////////////////////////////////////////////////

// STEPS library.
using steps::rng::R123;

/// Increment the first 2 values of ctr, translated to 64-bit
void ctr_increment(R123::r123_type::ctr_type &ctr) {
    uint64_t x = ctr[0] + (static_cast<uint64_t>(ctr[1]) << 32);
    x++;
    ctr[0] = x;
    ctr[1] = x >> 32;
}

////////////////////////////////////////////////////////////////////////////////
void R123::concreteInitialize(ulong seed)
{
    key.fill(0);
    ctr[0] = 0;          /// Incrementing the counter
    ctr[1] = 0;          /// Incrementing the counter
    ctr[2] = seed;       /// First 32 bits of 64-bit seed
    ctr[3] = seed >> 32; /// Last 32 bits of the seed
}

////////////////////////////////////////////////////////////////////////////////

/// Fills the buffer with random numbers on [0,0xffffffff]-interval.
void R123::concreteFillBuffer()
{
    uint* b;
    for (b = rBuffer; b+4 <= rEnd; b+=4) {
        /// Getting 4 new random numbers
        r123_type::ctr_type rn = r(ctr, key);
        ctr_increment(ctr);
        b[0] = rn[0];
        b[1] = rn[1];
        b[2] = rn[2];
        b[3] = rn[3];
    }

    if (b == rEnd) {
        return;
    }

    AssertLog(b+4 > rEnd);
    r123_type::ctr_type rn = r(ctr, key);
    ctr_increment(ctr);
    for (int i = 0; b < rEnd; ++b, ++i) {
        *b = rn[i];
    }
}

////////////////////////////////////////////////////////////////////////////////

// END
