/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2020 Okinawa Institute of Science and Technology, Japan.
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


#ifndef STEPS_RNG_RNG_HPP
#define STEPS_RNG_RNG_HPP 1


// STL headers.
#include <memory>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/tools.hpp"

namespace steps {
namespace rng {

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
    virtual ~RNG();

    /// Initialize the generator with seed.
    ///
    /// \param seed Seed for the generator.
    void initialize(ulong const & seed);

    /// Minimax inclusive range for the C++11 compatibility
    static constexpr uint min() { return 0; }
    static constexpr uint max() { return 0xffffffffu; }

    typedef uint result_type;
    result_type operator()() { return get(); }

    /// Return the next random int in the buffer of the generator.
    ///
    inline uint get()
    {
        if (rNext == rEnd) { concreteFillBuffer(); rNext = rBuffer; }
        return *(rNext++);
    }

    /// Generates a uniform random number on [0,1] real interval.
    ///
    inline double getUnfII()
    {
        // Divided by 2^32-1.
        return get() * (1.0 / 4294967295.0);
    }

    /// Generates a uniform random number on [0,1) real interval.
    ///
    inline double getUnfIE()
    {
        // Divided by 2^32.
        return get() * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on (0,1) real interval.
    ///
    inline double getUnfEE()
    {
        // Divided by 2^32.
        return (static_cast<double>(get()) + 0.5) * (1.0/4294967296.0);
    }

    /// Generates a uniform random number on [0,1) with 53-bit resolution.
    ///
    inline double getUnfIE53()
    {
        ulong a = get() >> 5, b = get() >> 6;
        return(a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
    }

    /// Get a standard exponentially distributed number.
    float getStdExp();

    /// Get an exponentially distributed number with mean lambda.
    ///
    double getExp(double lambda);

    /// Get a Poisson-distributed number with mean lambda.
    ///
    long getPsn(float lambda);

    /// Get a standard normally distributed random number.
    ///
    float getStdNrm();

    /// Get a binomially distributed number with parameters t and p.
    ///
    uint getBinom(uint t, double p);

protected:

    uint                      * rBuffer;
    uint                        rSize;

    uint                      * rNext;
    uint                      * rEnd;

    virtual void concreteInitialize(ulong seed) = 0;

    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
    ///
    virtual void concreteFillBuffer() = 0;

private:

    bool                        pInitialized;

};

using RNGptr = std::shared_ptr<RNG>;

////////////////////////////////////////////////////////////////////////////////

} // namespace rng
} // namespace steps

#endif
// STEPS_RNG_RNG_HPP

// END

