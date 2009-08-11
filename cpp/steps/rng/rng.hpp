////////////////////////////////////////////////////////////////////////////////
// STEPS - STochastic Engine for Pathway Simulation
// Copyright (C) 2005-2008 Stefan Wils. All rights reserved.
//
// This file is part of STEPS.
//
// This library is free software; you can redistribute it and/or
// modify it under the terms of the GNU Lesser General Public
// License as published by the Free Software Foundation; either
// version 2.1 of the License, or (at your option) any later version.
//
// This library is distributed in the hope that it will be useful,
// but WITHOUT ANY WARRANTY; without even the implied warranty of
// MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
// Lesser General Public License for more details.
//
// You should have received a copy of the GNU Lesser General Public
// License along with this library; if not, write to the Free Software
// Foundation, Inc., 51 Franklin St, Fifth Floor, Boston, MA  02110-1301, USA
//
//
////////////////////////////////////////////////////////////////////////////////


#ifndef STEPS_RNG_RNG_HPP
#define STEPS_RNG_RNG_HPP 1

// Autotools definitions.
#ifdef HAVE_CONFIG_H
#include <steps/config.h>
#endif

// STL headers.
#include <string>

// STEPS headers.
#include <steps/common.h>
#include <steps/math/tools.hpp>

START_NAMESPACE(steps)
START_NAMESPACE(rng)

////////////////////////////////////////////////////////////////////////////////

class RNG
{

public:

    RNG(uint bufsize);
    virtual ~RNG(void);

    void initialize(ulong const & seed);

    inline uint get(void)
    {
        if (rNext == rEnd) { concreteFillBuffer(); rNext = rBuffer; }
        return *(rNext++);
    }

    /// Generates a uniform random number on [0,1] real interval.
    inline double getUnfII(void)
    {
        // Divided by 2^32-1.
        return get() * (1.0 / 4294967295.0);
    }

    /// Generates a uniform random number on [0,1) real interval.
    inline double getUnfIE(void)
    {
        // Divided by 2^32.
        return get() * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on (0,1) real interval.
    inline double getUnfEE(void)
    {
        // Divided by 2^32.
        return (((double)get()) + 0.5) * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on [0,1) with 53-bit resolution.
    inline double getUnfIE53(void)
    {
        ulong a = get() >> 5, b = get() >> 6;
        return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

    /// Get a standard exponentially distributed number.
    float getStdExp(void);

    /// Get an exponentially distributed number with mean lambda.
    ///
    virtual double getExp(double lambda);

    /// Get a Poisson-distributed number with mean lambda.
    long getPsn(float lambda);

    /// Get a standard normally distributed random number.
    float getStdNrm(void);

protected:

    uint                      * rBuffer;
    uint                        rSize;

    uint                      * rNext;
    uint                      * rEnd;

    virtual void concreteInitialize(ulong seed) = 0;

    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
    virtual void concreteFillBuffer(void) = 0;

private:

    bool                        pInitialized;

};

////////////////////////////////////////////////////////////////////////////////

RNG * create_mt19937(uint bufsize);

RNG * create(std::string rng_name, uint bufsize);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(rng)
END_NAMESPACE(steps)

#endif
// STEPS_RNG_RNG_HPP

// END

