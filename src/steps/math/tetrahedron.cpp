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

#include "tetrahedron.hpp"

#include <cmath>
#include <limits>


#include "linsolve.hpp"
#include "point.hpp"
#include "triangle.hpp"

#include "util/error.hpp"

namespace steps::math {

////////////////////////////////////////////////////////////////////////////////

double tet_vol(const point3d& p0, const point3d& p1, const point3d& p2, const point3d& p3) {
    return std::abs((p1 - p0).dot((p2 - p0).cross(p3 - p0))) / 6;
}

////////////////////////////////////////////////////////////////////////////////

point3d tet_barycenter(const point3d& p0, const point3d& p1, const point3d& p2, const point3d& p3) {
    return (p0 + p1 + p2 + p3) / 4;
}

////////////////////////////////////////////////////////////////////////////////

// For a given face of a tetrahedron tell if a particular point (pi) is in the
// same side than the opposite point (opposite) than the triangle defined by
// (center, p1, p2). To know it, it only look if dot-product of (norm of
// (center, p1, p2) * (center, opposite) is the same sign than norm of (center,
// p1, p2) * (center, pi).
static inline bool same_direction(const point3d& center,
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
    auto tol = std::max(std::abs(A), std::abs(B)) * epsilon;

    // same sign with tolerance
    if (A < -tol && B > tol) {
        return false;
    }
    if (A > tol && B < -tol) {
        return false;
    }

    return true;
}

// TODO: try barycentic coordinates, it might be faster
bool tet_inside(const point3d& p0,
                const point3d& p1,
                const point3d& p2,
                const point3d& p3,
                const point3d& pi) {
    if (!same_direction(p0, p3, p1, p2, pi)) {
        return false;
    }
    if (!same_direction(p0, p2, p1, p3, pi)) {
        return false;
    }
    if (!same_direction(p0, p1, p2, p3, pi)) {
        return false;
    }
    if (!same_direction(p1, p0, p2, p3, pi)) {
        return false;
    }

    return true;
}

////////////////////////////////////////////////////////////////////////////////

point3d tet_ranpnt(const point3d& p0,
                   const point3d& p1,
                   const point3d& p2,
                   const point3d& p3,
                   double s,
                   double t,
                   double u) {
    // from
    // http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html
    if (s + t > 1.0) {
        // cut'n fold the cube into a prism
        s = 1.0 - s;
        t = 1.0 - t;
    }

    if (t + u > 1.0) {
        // cut'n fold the prism into a tetrahedron
        double tmp = u;
        u = 1.0 - s - t;
        t = 1.0 - tmp;
    } else if (s + t + u > 1.0) {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1.0 - t - tmp;
    }

    // a,s,t,u are the barycentric coordinates of the random point.
    double a = 1.0 - s - t - u;
    return a * p0 + s * p1 + t * p2 + u * p3;
}

////////////////////////////////////////////////////////////////////////////////

double tet_circumrad(const point3d& p0, const point3d& p1, const point3d& p2, const point3d& p3) {
    return std::sqrt(tet_circumrad2(p0, p1, p2, p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_circumrad2(const point3d& p0, const point3d& p1, const point3d& p2, const point3d& p3) {
    // Translate to p0 at origin
    point3d q1 = p1 - p0;
    point3d q2 = p2 - p0;
    point3d q3 = p3 - p0;

    double a[12];

    a[0] = q1[0];
    a[3] = q1[1];
    a[6] = q1[2];
    a[9] = q1.dot(q1);

    a[1] = q2[0];
    a[4] = q2[1];
    a[7] = q2[2];
    a[10] = q2.dot(q2);

    a[2] = q3[0];
    a[5] = q3[1];
    a[8] = q3[2];
    a[11] = q3.dot(q3);

    int info = math::linsolve(3, 1, a);
    if (info != 0) {
        AssertLog(0);
    }

    point3d d{a[9], a[10], a[11]};
    return 0.25 * d.dot(d);
}

////////////////////////////////////////////////////////////////////////////////

double tet_shortestedge(const point3d& p0,
                        const point3d& p1,
                        const point3d& p2,
                        const point3d& p3) {
    return std::sqrt(tet_shortestedge2(p0, p1, p2, p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_shortestedge2(const point3d& p0,
                         const point3d& p1,
                         const point3d& p2,
                         const point3d& p3) {
    return std::min(std::min(p0.dist2(p1), std::min(p0.dist2(p2), p0.dist2(p3))),
                    std::min(p1.dist2(p2), std::min(p1.dist2(p3), p2.dist2(p3))));
}

////////////////////////////////////////////////////////////////////////////////

double tet_longestedge(const point3d& p0, const point3d& p1, const point3d& p2, const point3d& p3) {
    return std::sqrt(tet_longestedge2(p0, p1, p2, p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_longestedge2(const point3d& p0,
                        const point3d& p1,
                        const point3d& p2,
                        const point3d& p3) {
    return std::max(std::max(p0.dist2(p1), std::max(p0.dist2(p2), p0.dist2(p3))),
                    std::max(p1.dist2(p2), std::max(p1.dist2(p3), p2.dist2(p3))));
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace steps::math
