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

#include <iosfwd>
#include <memory>

#include "math/tools.hpp"

namespace steps::rng {

////////////////////////////////////////////////////////////////////////////////
/// Base class of random number generator.
///
/// The RNG class can be inherited by other classes of random number generators.
class RNG {
  public:
    /// Constructor
    ///
    /// \param bufsize Size of the buffer.
    RNG(unsigned int bufsize);

    /// Destructor
    virtual ~RNG();

    virtual void checkpoint(std::ostream& cp_file) const;

    virtual void restore(std::istream& cp_file);

    /// Initialize the generator with seed.
    ///
    /// \param seed Seed for the generator.
    void initialize(unsigned long const& seed);

    /// Return the seed of the generator.
    unsigned long seed() const;

    /// Minimax inclusive range for the C++11 compatibility
    static constexpr unsigned int min() {
        return 0;
    }
    static constexpr unsigned int max() {
        return 0xffffffffu;
    }

    typedef unsigned int result_type;
    result_type operator()() {
        return get();
    }

    /// Return the next random int in the buffer of the generator.
    ///
    inline unsigned int get() {
        if (rNext == rEnd) {
            concreteFillBuffer();
            rNext = rBuffer.get();
        }
        return *(rNext++);
    }

    /// Generates a uniform random number on [0,1] real interval.
    ///
    inline double getUnfII() {
        // Divided by 2^32-1.
        return get() * (1.0 / 4294967295.0);
    }

    /// Generates a uniform random number on [0,1) real interval.
    ///
    inline double getUnfIE() {
        // Divided by 2^32.
        return get() * (1.0 / 4294967296.0);
    }

    /// Generates a uniform random number on (0,1) real interval.
    ///
    inline double getUnfEE() {
        // Divided by 2^32.
        return (static_cast<double>(get()) + 0.5) * (1.0 / 4294967296.0);
    }

    /// Generates a uniform random number on [0,1) with 53-bit resolution.
    ///
    inline double getUnfIE53() {
        unsigned long a = get() >> 5, b = get() >> 6;
        return (a * 67108864.0 + b) * (1.0 / 9007199254740992.0);
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
    unsigned int getBinom(unsigned int t, double p);

    /**
     * Function that rounds a double value to one of the closest straddling integers
     * with a probability dependent on the proximity
     * \tparam T the type of the integer returned
     * \param value to round to one of the closest integer
     */
    template <typename T>
    T stochastic_round(double value) {
        std::uniform_real_distribution<double> uniform(0.0, 1.0);
        const auto floored_value = std::floor(value);
        return static_cast<T>(floored_value) + ((value - floored_value) > uniform(*this) ? 1 : 0);
    }

    /**
     * Function that rounds a double value to one of the closest straddling integers
     * with a probability dependent on the proximity \tparam T the type of the
     * integer returned
     * \param value to round to one of the closest integer
     * \param upper_bound is the maximum value allowed for the rounded value
     */

    template <typename T>
    T stochastic_round(double value, T upper_bound) {
        return std::min(upper_bound, this->stochastic_round<T>(value));
    }

  protected:
    std::unique_ptr<unsigned int[]> rBuffer;
    unsigned int rSize;

    unsigned int* rNext;
    unsigned int* rEnd;

    virtual void concreteInitialize(unsigned long seed) = 0;

    /// Fills the buffer with random numbers on [0,0xffffffff]-interval.
    ///
    virtual void concreteFillBuffer() = 0;

  private:
    bool pInitialized;
    unsigned long pSeed{};
};

using RNGptr = std::shared_ptr<RNG>;

}  // namespace steps::rng
