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

#ifndef STEPS_MATH_POINT_HPP
#define STEPS_MATH_POINT_HPP 1

#include <array>
#include <cmath>

namespace steps {
namespace math {

/** Array of 3 double values representing 3-d point.
 *
 * (consider replacing / augmenting with SIMD-friendly
 * representation.)
 */

struct point3d: public std::array<double,3> {
    point3d() = default;

    explicit point3d(double x0) {
        data()[0]=x0;
    }

    point3d(double x0, double x1) {
        data()[0]=x0;
        data()[1]=x1;
    }

    point3d(double x0, double x1, double x2) {
        data()[0]=x0;
        data()[1]=x1;
        data()[2]=x2;
    }

    // arithmetic operations

    point3d &operator+=(const point3d &p) {
        (*this)[0]+=p[0];
        (*this)[1]+=p[1];
        (*this)[2]+=p[2];
        return *this;
    }

    point3d &operator-=(const point3d &p) {
        (*this)[0]-=p[0];
        (*this)[1]-=p[1];
        (*this)[2]-=p[2];
        return *this;
    }

    point3d &operator*=(double x) {
        (*this)[0]*=x;
        (*this)[1]*=x;
        (*this)[2]*=x;
        return *this;
    }

    point3d &operator/=(double x) {
        (*this)[0]/=x;
        (*this)[1]/=x;
        (*this)[2]/=x;
        return *this;
    }

    // equality testing

    bool operator==(const point3d &x) const {
        return (*this)[0]==x[0] && (*this)[1]==x[1] && (*this)[2]==x[2];
    }

    bool almostEqual(const point3d &x) const {
        static const double EPSILON = 1e-9;
        if (*this == x){
            return true;
        }
        const point3d d(std::abs((*this)[0] - x[0]),
                        std::abs((*this)[1] - x[1]),
                        std::abs((*this)[2] - x[2]));
        return d[0] < std::abs(x[0]) * EPSILON &&
               d[1] < std::abs(x[1]) * EPSILON &&
               d[2] < std::abs(x[2]) * EPSILON;
    }

    bool operator!=(const point3d &x) const {
        return !(*this==x);
    }

    static point3d zero() { return point3d{0,0,0}; }

    inline double dot(const point3d &q) const {
        const point3d &p = *this;
        return p[0]*q[0]+p[1]*q[1]+p[2]*q[2];
    }

    inline point3d cross(const point3d &q) const {
        const point3d &p = *this;
        return point3d{
                p[1]*q[2] - p[2]*q[1],
                p[2]*q[0] - p[0]*q[2],
                p[0]*q[1] - p[1]*q[0]
        };
    }

    inline double dist2(point3d b) const {
        b -= (*this);
        return b.dot(b);
    }

    inline double distance(const point3d &b) const {
        return std::sqrt(this->dist2(b));
    }

};

inline point3d operator+(point3d p, const point3d &q) {
    return p+=q;
}

inline point3d operator-(point3d p, const point3d &q) {
    return p-=q;
}

inline point3d operator*(point3d p, double x) {
    return p*=x;
}

inline point3d operator*(double x, point3d p) {
    return p*=x;
}

inline point3d operator/(point3d p, double x) {
    return p/=x;
}


// Compat. point3d overload
inline double distance(const point3d &a, const point3d &b) {
    return a.distance(b);
}

}} // namespace steps::math

#endif // ndef STEPS_MATH_POINT_HPP
