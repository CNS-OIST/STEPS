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

#include "triangle.hpp"
#include "segment.hpp"

#include <cassert>
#include <cmath>
#include <limits>

namespace steps::math {

////////////////////////////////////////////////////////////////////////////////

double tri_area(const point3d& p0, const point3d& p1, const point3d& p2) {
    point3d c = (p1 - p0).cross(p2 - p0);
    return 0.5 * norm(c);
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_barycenter(const point3d& p0, const point3d& p1, const point3d& p2) {
    return (p0 + p1 + p2) / 3.0;
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_normal(const point3d& p0, const point3d& p1, const point3d& p2) {
    point3d edge1 = p1 - p0;
    point3d edge2 = p2 - p0;
    point3d c = edge1.cross(edge2);
    const auto twice_area = norm(c);
    assert(twice_area > 0);
    return c / twice_area;
}

////////////////////////////////////////////////////////////////////////////////

point3d tri_ranpnt(const point3d& p0, const point3d& p1, const point3d& p2, double s, double t) {
    // from http://www.cs.princeton.edu/~funk/tog02.pdf
    double u = std::sqrt(s);
    double v = u * t;

    return (1.0 - u) * p0 + (u - v) * p1 + v * p2;
}

////////////////////////////////////////////////////////////////////////////////

bool tri_intersect_line(const point3d& v0,
                        const point3d& v1,
                        const point3d& v2,
                        const point3d& p0,
                        const point3d& p1,
                        point3d& intersection,
                        bool is_segment) {
    // Compute the intersection point using moller-trumbore method
    // https://en.wikipedia.org/wiki/M%C3%B6ller%E2%80%93Trumbore_intersection_algorithm
    //
    // Notation follows mostly the paper:
    //
    // Moeller, Tomas; Trumbore, Ben (1997).
    // "Fast, Minimum Storage Ray-Triangle Intersection".
    // Journal of Graphics Tools. 2: 21–28
    //
    // u,v are the triangle baricentric coordinates
    // t is the parametric distance along the ray

    const auto edge1 = v1 - v0, edge2 = v2 - v0;
    const auto pdir = p1 - p0;

    const auto pvec = pdir.cross(edge2);
    const auto det = edge1.dot(pvec);

    // 4 ULP
    constexpr auto tol = 4 * std::numeric_limits<point3d::value_type>::epsilon();
    // get max volume to scale tolerance (det = scalar triple product)
    // use cheaper L1 norm since |L1| >= |L2|
    const auto max_volume = normL1(pdir) * normL1(edge1) * normL1(edge2);
    const auto tol_volume = tol * max_volume;

    // line and triangle parallel?
    if (det > -tol_volume && det < tol_volume) {
        return false;
    }

    const auto invDet = 1.0 / det;

    // calculate distance from v0 and ray origin
    auto tvec = p0 - v0;

    // calculate u coo and test bounds
    const auto u = tvec.dot(pvec) * invDet;
    if (u < -tol || u > 1.0 + tol) {
        return false;
    }

    // calculate v coo and test bounds
    const auto qvec = tvec.cross(edge1);
    const auto v = pdir.dot(qvec) * invDet;
    if (v < -tol || u + v > 1.0 + tol) {
        return false;
    }

    // calculate t coo and test bounds
    const auto t = edge2.dot(qvec) * invDet;
    if (t < -tol || (is_segment && t > 1.0 + tol)) {
        return false;
    }
    intersection = p0 + t * pdir;
    return true;
}


////////////////////////////////////////////////////////////////////////////////


bool tri_intersect_point(const point3d& v0,
                         const point3d& v1,
                         const point3d& v2,
                         const point3d& p,
                         const bool is_triangle) {
    // barycentric coordinates method
    const auto edge0 = v1 - v0;
    const auto edge1 = v2 - v0;

    // degenerate triangles
    if (segment_intersect_point(v0, v1, v2, false)) {
        // it is ok to use the L1 norm because the triangle is colinear
        const auto n0 = normL1(edge0);
        const auto n1 = normL1(edge1);
        const auto n2 = normL1(v2 - v1);
        if (n0 >= n1 && n0 >= n2) {
            return segment_intersect_point(v0, v1, p);
        } else if (n1 >= n0 && n1 >= n2) {
            return segment_intersect_point(v2, v0, p);
        } else if (n2 >= n0 && n2 >= n1) {
            return segment_intersect_point(v2, v1, p);
        }
    }

    // 4ULP
    constexpr auto adimensional_tol = 4 * std::numeric_limits<point3d::value_type>::epsilon();
    const auto tol = adimensional_tol *
                     std::max(std::max(normL1(edge0), normL1(edge1)), normL1(v2 - v1));

    if (!is_triangle) {
        // just test if point p belongs to surface defined by the triangle
        const auto normal = (edge0).cross(edge1);
        const auto A = std::abs((p - v0).dot(normal));
        return A < tol;
    }

    const auto pvec = p - v0;

    const double d00 = edge0.dot(edge0);
    const double d01 = edge0.dot(edge1);
    const double d11 = edge1.dot(edge1);
    const double d20 = pvec.dot(edge0);
    const double d21 = pvec.dot(edge1);

    const double denom = d00 * d11 - d01 * d01;

    // barycentric coordinates multiplied by denom to avoid divisions
    const double c1 = (d11 * d20 - d01 * d21);
    const double c2 = (d00 * d21 - d01 * d20);
    const double c0 = denom - c1 - c2;

    const double denomtol = denom * adimensional_tol;

    // check that the p projection of the plane is inside the triangle
    if (c0 < -denomtol || c0 > denom + denomtol) {
        return false;
    }
    if (c1 < -denomtol || c1 > denom + denomtol) {
        return false;
    }
    if (c2 < -denomtol || c2 > denom + denomtol) {
        return false;
    }

    // check that projection and initial point overlap
    const auto pproj = (c0 / denom) * v0 + (c1 / denom) * v1 + (c2 / denom) * v2;
    return p.almostEqual(pproj, tol);
}

////////////////////////////////////////////////////////////////////////////////

bool same_direction(const point3d& center,
                    const point3d& opposite,
                    const point3d& p1,
                    const point3d& p2,
                    const point3d& pi) {
    // TetMesh check that any tetrahedron is not degenerated, so don't test it

    // FIXME: the normal of a face is alreay computed by TetMesh, but I don't see
    // how to use it.
    auto normal = (p1 - center).cross(p2 - center);
    auto A = (opposite - center).dot(normal);
    auto B = (pi - center).dot(normal);

    // 4 ULP
    constexpr auto epsilon = 4 * std::numeric_limits<point3d::value_type>::epsilon();
    const double tol = std::max(std::abs(A), std::abs(B)) * epsilon;

    // same sign with tolerance
    if (A < -tol && B > tol) {
        return false;
    }
    if (A > tol && B < -tol) {
        return false;
    }

    return true;
}


}  // namespace steps::math
