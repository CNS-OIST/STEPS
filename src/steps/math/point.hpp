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

#include <array>
#include <cmath>
#include <limits>

namespace steps::math {

// maximum linear geometric scale expected
// ideally this should be set as max(abs(x,y,z)) from the mesh
static constexpr double max_length_scale = 1000.0;
// linear tolerance (4 ULP)
static constexpr double tol_lin = 4 * max_length_scale * std::numeric_limits<double>::epsilon();

struct point3d_trait {};
struct position_abs_trait {};
struct position_rel_to_ves_trait {};

template <typename TFrom, typename TTo>
struct is_conversion: std::false_type {};

template <>
struct is_conversion<point3d_trait, position_abs_trait>: std::true_type {};
template <>
struct is_conversion<position_abs_trait, point3d_trait>: std::true_type {};

/** Array of 3 double values representing 3-d point.
 *
 * (consider replacing / augmenting with SIMD-friendly
 * representation.)
 */

template <typename Trait>
struct point3d_: public std::array<double, 3> {
    point3d_() {
        data()[0] = 0.0;
        data()[1] = 0.0;
        data()[2] = 0.0;
    }

    point3d_(double x0, double x1, double x2) {
        data()[0] = x0;
        data()[1] = x1;
        data()[2] = x2;
    }

    template <typename OTrait,
              typename = typename std::enable_if<is_conversion<OTrait, Trait>::value>::type>
    point3d_(const point3d_<OTrait>& other) {
        data()[0] = other[0];
        data()[1] = other[1];
        data()[2] = other[2];
    }

    // arithmetic operations

    point3d_<Trait>& operator+=(const point3d_<Trait>& p) {
        (*this)[0] += p[0];
        (*this)[1] += p[1];
        (*this)[2] += p[2];
        return *this;
    }

    point3d_<Trait>& operator-=(const point3d_<Trait>& p) {
        (*this)[0] -= p[0];
        (*this)[1] -= p[1];
        (*this)[2] -= p[2];
        return *this;
    }

    point3d_<Trait>& operator*=(double x) {
        (*this)[0] *= x;
        (*this)[1] *= x;
        (*this)[2] *= x;
        return *this;
    }

    point3d_<Trait>& operator/=(double x) {
        (*this)[0] /= x;
        (*this)[1] /= x;
        (*this)[2] /= x;
        return *this;
    }

    // equality testing

    bool operator==(const point3d_<Trait>& x) const {
        return (*this)[0] == x[0] && (*this)[1] == x[1] && (*this)[2] == x[2];
    }

    bool almostEqual(const point3d_<Trait>& x) const {
        if (*this == x) {
            return true;
        }
        const point3d_<Trait> d(std::abs((*this)[0] - x[0]),
                                std::abs((*this)[1] - x[1]),
                                std::abs((*this)[2] - x[2]));
        return d[0] < tol_lin && d[1] < tol_lin && d[2] < tol_lin;
    }

    bool operator!=(const point3d_<Trait>& x) const {
        return !(*this == x);
    }

    static constexpr point3d_<Trait> zero() {
        return point3d_<Trait>{0, 0, 0};
    }

    inline double dot(const point3d_<Trait>& q) const {
        const point3d_<Trait>& p = *this;
        return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
    }

    inline point3d_<Trait> cross(const point3d_<Trait>& q) const {
        const point3d_<Trait>& p = *this;
        return point3d_<Trait>{p[1] * q[2] - p[2] * q[1],
                               p[2] * q[0] - p[0] * q[2],
                               p[0] * q[1] - p[1] * q[0]};
    }

    inline double dist2(point3d_<Trait> b) const {
        b -= (*this);
        return b.dot(b);
    }

    inline double distance(const point3d_<Trait>& b) const {
        return std::sqrt(this->dist2(b));
    }
};

template <typename Trait>
inline point3d_<Trait> operator+(point3d_<Trait> p, const point3d_<Trait>& q) {
    return p += q;
}

template <typename Trait>
inline point3d_<Trait> operator-(point3d_<Trait> p, const point3d_<Trait>& q) {
    return p -= q;
}

template <typename Trait>
inline point3d_<Trait> operator*(point3d_<Trait> p, double x) {
    return p *= x;
}

template <typename Trait>
inline point3d_<Trait> operator*(double x, point3d_<Trait> p) {
    return p *= x;
}

template <typename Trait>
inline point3d_<Trait> operator/(point3d_<Trait> p, double x) {
    return p /= x;
}

template <typename Trait>
inline double mag(const point3d_<Trait>& p) {
    return std::sqrt(p[0] * p[0] + p[1] * p[1] + p[2] * p[2]);
}

template <typename Trait>
inline double dot(const point3d_<Trait>& p, const point3d_<Trait>& q) {
    return p[0] * q[0] + p[1] * q[1] + p[2] * q[2];
}

template <typename Trait>
inline point3d_<Trait> cross(const point3d_<Trait>& p, const point3d_<Trait>& q) {
    return point3d_<Trait>{p[1] * q[2] - p[2] * q[1],
                           p[2] * q[0] - p[0] * q[2],
                           p[0] * q[1] - p[1] * q[0]};
}

template <typename Trait>
inline double dist2(point3d_<Trait> a, const point3d_<Trait>& b) {
    a -= b;
    return dot(a, a);
}

// Compat. point3d overload
template <typename Trait>
inline double distance(const point3d_<Trait>& a, const point3d_<Trait>& b) {
    return a.distance(b);
}

template <typename Trait>
inline double squared_norm(const point3d_<Trait>& p) {
    return p.dot(p);
}

// L2-norm
template <typename Trait>
inline double norm(const point3d_<Trait>& p) {
    return std::sqrt(squared_norm(p));
}

// L1-norm, aka manhattan distance
template <typename Trait>
inline double normL1(const point3d_<Trait>& p) {
    return std::abs(p[0]) + std::abs(p[1]) + std::abs(p[2]);
}

// 3D point
using point3d = point3d_<point3d_trait>;
// An absolute position in space
using position_abs = point3d_<position_abs_trait>;
// A relative position to vesicle centre
using position_rel_to_ves = point3d_<position_rel_to_ves_trait>;

// Spherical coordinates
struct position_spherical: public std::array<double, 3> {
    position_spherical() {
        data()[0] = 0.0;
        data()[1] = 0.0;
        data()[2] = 0.0;
    }

    position_spherical(double radius, double theta, double phi) {
        data()[0] = radius;
        data()[1] = theta;
        data()[2] = phi;
    }

    double getRadius() {
        return data()[0];
    }
    double getTheta() {
        return data()[1];
    }
    double getPhi() {
        return data()[2];
    }

    void setRadius(double radius) {
        data()[0] = radius;
    }
    void setTheta(double theta) {
        data()[1] = theta;
    }
    void setPhi(double phi) {
        data()[2] = phi;
    }
};

}  // namespace steps::math
