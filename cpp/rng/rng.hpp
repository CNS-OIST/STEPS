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

#ifndef STEPS_RNG_RNG_HPP
#define STEPS_RNG_RNG_HPP 1


// STL headers.
#include <string>

// STEPS headers.
#include "../common.h"
#include "../math/tools.hpp"

START_NAMESPACE(steps)
START_NAMESPACE(rng)

////////////////////////////////////////////////////////////////////////////////
/// Base class of random number generator.
///
/// The RNG class can be inherited by other classes of random number generators.
class RNG
{

public:

    /// Constructor
    ///
    /// \param bufsize Size of the buffer.
    RNG(uint bufsize);

    /// Destructor
    virtual ~RNG(void);

    /// Initialize the generator with seed.
    ///
    /// \param seed Seed for the generator.
    void initialize(ulong const & seed);

    /// Return the next random int in the buffer of the generator.
    ///
    inline uint get(void)
    {
        if (rNext == rEnd) { concreteFillBuffer(); rNext = rBuffer; }
        return *(rNext++);
    }

    /// Generates a uniform random number on [0,1] real interval.
    ///
    inline double getUnfII(void)
    {
        // Divided by 2^32-1.
        return get() * (1.0 / 4294967295.0);
    }

    /// Generates a uniform random number on [0,1) real interval.
    ///
    inline double getUnfIE(void)
    {
        // Divided by 2^32.
        return get() * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on (0,1) real interval.
    ///
    inline double getUnfEE(void)
    {
        // Divided by 2^32.
        return (((double)get()) + 0.5) * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on [0,1) with 53-bit resolution.
    ///
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
    ///
    long getPsn(float lambda);

    /// Get a standard normally distributed random number.
    ///
    float getStdNrm(void);

protected:

    uint                      * rBuffer;
    uint                        rSize;

    uint                      * rNext;
    uint                      * rEnd;

    virtual void concreteInitialize(ulong seed) = 0;

    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
    ///
    virtual void concreteFillBuffer(void) = 0;

private:

    bool                        pInitialized;

};

////////////////////////////////////////////////////////////////////////////////
/// Create a MT19937 random number generator and return as RNG object.
///
/// \param buffsize Size of buffer.
RNG * create_mt19937(uint bufsize);

/// Create a random number generator with name rng_name and return as RNG object.
///
/// \param rng_name Name of the random number generator.
/// \param buffsize Size of buffer.

RNG * create(std::string rng_name, uint bufsize);

////////////////////////////////////////////////////////////////////////////////

END_NAMESPACE(rng)
END_NAMESPACE(steps)

#endif
// STEPS_RNG_RNG_HPP

// END

