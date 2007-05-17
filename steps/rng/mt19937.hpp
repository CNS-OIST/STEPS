////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2007 Stefan Wils. All rights reserved.
//
// $Id$
////////////////////////////////////////////////////////////////////////////////

#ifndef STEPS_RNG_MT19937_HPP
#define STEPS_RNG_MT19937_HPP 1

// STEPS headers.
#include <steps/common.h>
#include <steps/rng/rng.hpp>

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

class MT19937
: public RNG
{

public:

    MT19937(uint bufsize);
    virtual ~MT19937(void);

protected:

    virtual void concreteInitialize(ulong seed);
    
    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
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
