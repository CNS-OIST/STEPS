////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2007-2009ÊOkinawa Institute of Science and Technology, Japan.
// Copyright (C) 2003-2006ÊUniversity of Antwerp, Belgium.
//
// See the file AUTHORS for details.
//
// This file is part of STEPS.
//
// STEPSÊis free software: you can redistribute it and/or modify
// it under the terms of the GNU General Public License as published by
// the Free Software Foundation, either version 3 of the License, or
// (at your option) any later version.
//
// STEPSÊis distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
// GNU General Public License for more details.
//
// You should have received a copy of the GNU General Public License
// along with this program.ÊIf not, see <http://www.gnu.org/licenses/>.
//
////////////////////////////////////////////////////////////////////////////////

/*
 *  Last Changed Rev:  $Rev: 307 $
 *  Last Changed Date: $Date: 2010-03-24 09:25:37 +0900 (Wed, 24 Mar 2010) $
 *  Last Changed By:   $Author: wchen $
 */

#ifndef STEPS_RNG_MT19937_HPP
#define STEPS_RNG_MT19937_HPP 1


// STEPS headers.
#include "../common.h"
#include "rng.hpp"

START_NAMESPACE(steps)
START_NAMESPACE(rng)

////////////////////////////////////////////////////////////////////////////////

// Based on the original Mersenne Twister code (mt19937.c).
//
// Copyright (C) 1997 - 2002, Makoto Matsumoto and Takuji Nishimura,
// All rights reserved.
//
// Redistribution and use in source and binary forms, with or without
// modification, are permitted provided that the following conditions
// are met:
//
//   1. Redistributions of source code must retain the above copyright
//      notice, this list of conditions and the following disclaimer.
//
//   2. Redistributions in binary form must reproduce the above copyright
//      notice, this list of conditions and the following disclaimer in the
//      documentation and/or other materials provided with the distribution.
//
//   3. The names of its contributors may not be used to endorse or promote
//      products derived from this software without specific prior written
//      permission.
//
// THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
// "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
// LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
// A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL THE COPYRIGHT OWNER OR
// CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
// EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
// PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
// PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
// LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
// NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
// SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//
// Any feedback is very welcome.
// http://www.math.sci.hiroshima-u.ac.jp/~m-mat/MT/emt.html
// email: m-mat @ math.sci.hiroshima-u.ac.jp (remove space)

////////////////////////////////////////////////////////////////////////////////

#define MT_N                                    624

#define MT_M                                    397

#define MT_MATRIX_A                             0x9908b0dfUL

#define MT_UPPER_MASK                           0x80000000UL

#define MT_LOWER_MASK                           0x7fffffffUL

////////////////////////////////////////////////////////////////////////////////

/// MT19937 Random number generator.
///
/// This random number generator is based on the original Mersenne Twister code (mt19937.c).

class MT19937
: public RNG
{

public:

    /// Constructor
    ///
    /// \param bufsize Size of the buffer.
    MT19937(uint bufsize);

    /// Destructor
    ///
    virtual ~MT19937(void);

protected:

    /// Initialize the generator with seed.
    ///
    /// \param seed Seed for the generator.
    virtual void concreteInitialize(ulong seed);

    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
    ///
    virtual void concreteFillBuffer(void);

private:

    unsigned long               pState[MT_N];
    int                         pStateInit;

};

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(rng)
END_NAMESPACE(steps)

#endif
// STEPS_RNG_MT19937_HPP

// END
