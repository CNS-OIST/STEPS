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

#include <cmath>
#include <cstdlib>
#include <random>
#include <vector>

#include "constants.hpp"
#include "rng/rng.hpp"

namespace steps::math {

////////////////////////////////////////////////////////////////////////////////

inline bool isLargerEps(float r) {
    if (std::abs(r) > static_cast<float>(IEEE_EPSILON32))
        return true;
    return false;
}

inline bool isLargerEps(double r) {
    if (std::abs(r) > static_cast<double>(IEEE_EPSILON64))
        return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

inline bool isSmallerEps(float r) {
    if (std::abs(r) <= static_cast<float>(IEEE_EPSILON32))
        return true;
    return false;
}

inline bool isSmallerEps(double r) {
    if (std::abs(r) <= static_cast<double>(IEEE_EPSILON64))
        return true;
    return false;
}

////////////////////////////////////////////////////////////////////////////////

inline float getSysRand(float min, float max) {
    return min + ((max - min) * static_cast<float>(rand()) / static_cast<float>(RAND_MAX));
}

inline double getSysRand(double min, double max) {
    return min + ((max - min) * static_cast<double>(rand()) / static_cast<double>(RAND_MAX));
}

inline int getSysRand(int min, int max) {
    return static_cast<int>(
        floorf(0.5 + getSysRand(static_cast<double>(min), static_cast<double>(max))));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T min(T const& v1, T const& v2) {
    return v1 < v2 ? v1 : v2;
}

template <typename T>
inline T min(T const& v1, T const& v2, T const& v3) {
    return min(v1, min(v2, v3));
}

template <typename T>
inline T min(T const& v1, T const& v2, T const& v3, T const& v4) {
    return min(min(v1, v2), min(v3, v4));
}

////////////////////////////////////////////////////////////////////////////////

template <typename T>
inline T max(T const& v1, T const& v2) {
    return v1 > v2 ? v1 : v2;
}

template <typename T>
inline T max(T const& v1, T const& v2, T const& v3) {
    return max(v1, max(v2, v3));
}

template <typename T>
inline T max(T const& v1, T const& v2, T const& v3, T const& v4) {
    return max(max(v1, v2), max(v3, v4));
}

////////////////////////////////////////////////////////////////////////////////

// Transfers sign of argument sign to argument num.
template <typename T>
inline T sign(T const& num, T const& sign) {
    if (((sign > 0) && (num < 0)) || ((sign < 0) && (num > 0)))
        return -num;
    else
        return num;
}

////////////////////////////////////////////////////////////////////////////////

extern void setSysRandInitTime();

////////////////////////////////////////////////////////////////////////////////

inline float xformDegToRad(float deg) {
    return (deg * static_cast<float>(PI)) / 180.0f;
}

inline double xformDegToRad(double deg) {
    return (deg * static_cast<double>(PI)) / 180.0;
}

////////////////////////////////////////////////////////////////////////////////

inline float xformRadToDeg(float rad) {
    return (rad * 180.0f) / static_cast<float>(PI);
}

inline double xformRadToDeg(double rad) {
    return (rad * 180.0) / static_cast<double>(PI);
}

////////////////////////////////////////////////////////////////////////////////

/// The factorial function
inline unsigned int factorial(unsigned int n) {
    if (n == 0)
        return 1;
    return n * factorial(n - 1);
}

/// The binomial coefficients
inline unsigned int binom_coeff(unsigned int n, unsigned int k) {
    //    return factorial(n)/factorial(k)/factorial(n-k);
    return round(exp(lgamma(n + 1) - lgamma(n - k + 1) - lgamma(k + 1)));
}

////////////////////////////////////////////////////////////////////////////////

inline void cross_product(double a1,
                          double a2,
                          double a3,
                          double b1,
                          double b2,
                          double b3,
                          double& p1,
                          double& p2,
                          double& p3) {
    p1 = (a2 * b3 - a3 * b2);
    p2 = (a3 * b1 - a1 * b3);
    p3 = (a1 * b2 - a2 * b1);
}

inline double dot_product(double a1, double a2, double a3, double b1, double b2, double b3) {
    return a1 * b1 + a2 * b2 + a3 * b3;
}

}  // namespace steps::math
