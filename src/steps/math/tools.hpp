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

#ifndef STEPS_MATH_TOOLS_HPP
#define STEPS_MATH_TOOLS_HPP 1


// Standard library & STL headers.
#include <cmath>
#include <cstdlib>
#include <vector>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/constants.hpp"

namespace steps {
namespace math {

////////////////////////////////////////////////////////////////////////////////

inline bool isLargerEps(float r)
{
    if (fabsf(r) > static_cast<float>(IEEE_EPSILON32)) return true;
    return false;
}

inline bool isLargerEps(double r)
{
    if (fabs(r) > static_cast<double>(IEEE_EPSILON64)) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

inline bool isSmallerEps(float r)
{
    if (fabsf(r) <= static_cast<float>(IEEE_EPSILON32)) return true;
    return false;
}

inline bool isSmallerEps(double r)
{
    if (fabs(r) <= static_cast<double>(IEEE_EPSILON64)) return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

inline float getSysRand(float min, float max)
{
    return min + ((max - min) *
        static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
}

inline double getSysRand(double min, double max)
{
    return min + ((max - min) *
        static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
}

inline int getSysRand(int min, int max)
{
    return static_cast<int>(floorf(0.5 + getSysRand(static_cast<double>(min),
                                              static_cast<double>(max))));
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
inline T min(T const & v1, T const & v2)
{
    return (v1 < v2 ? v1 : v2);
}

template<typename T>
inline T min(T const & v1, T const & v2, T const & v3)
{
    return min(v1, min(v2, v3));
}

template<typename T>
inline T min(T const & v1, T const & v2, T const & v3, T const & v4)
{
    return min(min(v1, v2), min(v3, v4));
}

////////////////////////////////////////////////////////////////////////////////

template<typename T>
inline T max(T const & v1, T const & v2)
{
    return (v1 > v2 ? v1 : v2);
}

template<typename T>
inline T max(T const & v1, T const & v2, T const & v3)
{
    return max(v1, max(v2, v3));
}

template<typename T>
inline T max(T const & v1, T const & v2, T const & v3, T const & v4)
{
    return max(max(v1, v2), max(v3, v4));
}

////////////////////////////////////////////////////////////////////////////////

// Transfers sign of argument sign to argument num.
template<typename T>
inline T sign(T const & num, T const & sign)
{
    if (((sign > 0) && (num < 0)) || ((sign < 0) && (num > 0)))
        return -num;
    else
        return num;
}

////////////////////////////////////////////////////////////////////////////////

extern void setSysRandInitTime();

////////////////////////////////////////////////////////////////////////////////

inline float xformDegToRad(float deg)
{
    return (deg * static_cast<float>(PI)) / 180.0f;
}

inline double xformDegToRad(double deg)
{
    return (deg * static_cast<double>(PI)) / 180.0;
}

////////////////////////////////////////////////////////////////////////////////

inline float xformRadToDeg(float rad)
{
    return (rad * 180.0f) / static_cast<float>(PI);
}

inline double xformRadToDeg(double rad)
{
    return (rad * 180.0) / static_cast<double>(PI);
}

////////////////////////////////////////////////////////////////////////////////

/// The factorial function
inline uint factorial(uint n)
{
    if (n == 0)
        return 1;
    return n*factorial(n-1);
}

/// The binomial coefficients
inline uint binom_coeff(uint n, uint k)
{
//    return factorial(n)/factorial(k)/factorial(n-k);
    return round(exp(lgamma(n+1) - lgamma(n-k+1) - lgamma(k+1)));
}


////////////////////////////////////////////////////////////////////////////////

}
}

#endif
// STEPS_MATH_TOOLS_HPP

// END
