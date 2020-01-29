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

#include <cmath>

// STEPS headers.
#include "steps/common.h"
#include "steps/math/point.hpp"
#include "steps/math/triangle.hpp"

namespace steps {
namespace math {

////////////////////////////////////////////////////////////////////////////////

double tri_area(const point3d &p0, const point3d &p1, const point3d &p2)
{
    point3d c = (p1-p0).cross(p2-p0);
    return 0.5 * std::sqrt(c.dot(c));
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_barycenter(const point3d &p0, const point3d &p1, const point3d &p2)
{
    return (p0+p1+p2)/3.0;
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_normal(const point3d &p0, const point3d &p1, const point3d &p2)
{
    point3d c = (p1-p0).cross(p2-p0);
    return c/std::sqrt(c.dot(c));
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_ranpnt(const point3d &p0, const point3d &p1, const point3d &p2, double s, double t)
{
    // from http://www.cs.princeton.edu/~funk/tog02.pdf
    double u = std::sqrt(s);
    double v = u*t;

    return (1-u)*p0 + (u-v)*p1 + v*p2;
}

////////////////////////////////////////////////////////////////////////////////

bool tri_intersect_line(const point3d &v0, const point3d &v1, const point3d &v2,
                        const point3d &p0, const point3d &p1, point3d &intersection,
                        bool is_segment)
{
    // Compute the intersection point using moller-trumbore method
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    using vec_t = point3d;
    const vec_t edge1 = v1 - v0,
                edge2 = v2 - v0;

    // EPSILON should be relative. We use a quick Manhattan Distance as reference
    const double EPSILON = std::abs(edge1[0] + edge1[1] + edge1[2]) * 1.e-5;

    const vec_t pdir = p1 - p0;
    const vec_t pvec = pdir.cross(edge2);
    const double det = edge1.dot(pvec);

    // line and triangle parallel?
    if (det > -EPSILON && det < EPSILON) return false;

    const double invDet = 1.0 / det;

    vec_t t = p0 - v0;
    double u = t.dot(pvec) * invDet;
    if (u < .0 || u > 1.0) return false;

    vec_t q = t.cross(edge1);
    double v = pdir.dot(q) * invDet;
    if (v < .0 || u + v > 1.0) return false;

    double x = edge2.dot(q) * invDet;
    if (x < EPSILON || (is_segment && x > 1.0)) return false;
    intersection = p0 + x * pdir;
    return true;
}

}  // namespace math
}  // namespace steps

// END
