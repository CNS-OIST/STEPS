/*
 #################################################################################
#
#    STEPS - STochastic Engine for Pathway Simulation
#    Copyright (C) 2007-2018 Okinawa Institute of Science and Technology, Japan.
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
#include "steps/error.hpp"
#include "steps/math/linsolve.hpp"
#include "steps/math/point.hpp"
#include "steps/math/tetrahedron.hpp"

// logging
#include "easylogging++.h"

namespace steps {
namespace math {

////////////////////////////////////////////////////////////////////////////////

static std::array<double,4> tet_barycentric(const point3d &p0, point3d p1,
                                            point3d p2, point3d p3, point3d pi) 
{
    // Translate to p0 at origin
    p1-=p0;
    p2-=p0;
    p3-=p0;
    pi-=p0;

    // Prepase base matrix.
    double a[12];
    a[ 0] = p1[0];
    a[ 1] = p1[1];
    a[ 2] = p1[2];
    a[ 3] = p2[0];
    a[ 4] = p2[1];
    a[ 5] = p2[2];
    a[ 6] = p3[0];
    a[ 7] = p3[1];
    a[ 8] = p3[2];
    a[ 9] = pi[0];
    a[10] = pi[1];
    a[11] = pi[2];

    // Solve this system.
    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
        AssertLog(0);

    double b1 = a[ 9];
    double b2 = a[10];
    double b3 = a[11];
    return std::array<double,4>{1 - (b1 + b2 + b3), b1, b2, b3};
}

////////////////////////////////////////////////////////////////////////////////

double tet_vol(const point3d &p0, const point3d &p1,
                            const point3d &p2, const point3d &p3)
{
    return std::fabs(dot(p1-p0,cross(p2-p0,p3-p0)))/6;
}

////////////////////////////////////////////////////////////////////////////////

point3d tet_barycenter(const point3d &p0, const point3d &p1,
                       const point3d &p2, const point3d &p3)
{
    return (p0+p1+p2+p3)/4;
}

////////////////////////////////////////////////////////////////////////////////

bool tet_inside(const point3d &p0, const point3d &p1,
                const point3d &p2, const point3d &p3,
                const point3d &pi)
{
    auto po = tet_barycentric(p0, p1, p2, p3, pi);
    return po[0] >= 0 && po[1] >= 0 && po[2] >= 0 && po[3] >= 0;
}

////////////////////////////////////////////////////////////////////////////////

point3d tet_ranpnt(const point3d &p0, const point3d &p1,
                   const point3d &p2, const point3d &p3,
                   double s, double t, double u)
{
    // from http://vcg.isti.cnr.it/activities/geometryegraphics/pointintetraedro.html
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
    }
    else if (s + t + u > 1.0) {
        double tmp = u;
        u = s + t + u - 1.0;
        s = 1 - t - tmp;
    }

    // a,s,t,u are the barycentric coordinates of the random point.
    double a = 1 - s - t - u;
    return a*p0 + s*p1 + t*p2 + u*p3;
}

////////////////////////////////////////////////////////////////////////////////

double tet_circumrad(const point3d &p0, const point3d &p1,
                     const point3d &p2, const point3d &p3)
{
    return std::sqrt(tet_circumrad2(p0,p1,p2,p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_circumrad2(const point3d &p0, const point3d &p1,
                      const point3d &p2, const point3d &p3)
{
    // Translate to p0 at origin
    point3d q1=p1-p0;
    point3d q2=p2-p0;
    point3d q3=p3-p0;

    double a[12];

    a[ 0] = q1[0];
    a[ 3] = q1[1];
    a[ 6] = q1[2];
    a[ 9] = dot(q1, q1);

    a[ 1] = q2[0];
    a[ 4] = q2[1];
    a[ 7] = q2[2];
    a[10] = dot(q2, q2);

    a[ 2] = q3[0];
    a[ 5] = q3[1];
    a[ 8] = q3[2];
    a[11] = dot(q3, q3);

    int info = steps::math::linsolve(3, 1, a);
    if (info != 0)
        AssertLog(0);

    point3d d{a[9], a[10], a[11]};
    return 0.25 * dot(d, d);
}

////////////////////////////////////////////////////////////////////////////////

double tet_shortestedge(const point3d &p0, const point3d &p1,
                        const point3d &p2, const point3d &p3)
{
    return std::sqrt(tet_shortestedge2(p0, p1, p2, p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_shortestedge2(const point3d &p0, const point3d &p1,
                         const point3d &p2, const point3d &p3)
{
    return std::min(std::min(dist2(p0,p1), std::min(dist2(p0,p2), dist2(p0,p3))),
                    std::min(dist2(p1,p2), std::min(dist2(p1,p3), dist2(p2,p3))));
}

////////////////////////////////////////////////////////////////////////////////

double tet_longestedge(const point3d &p0, const point3d &p1,
                       const point3d &p2, const point3d &p3)
{
    return std::sqrt(tet_longestedge2(p0, p1, p2, p3));
}

////////////////////////////////////////////////////////////////////////////////

double tet_longestedge2(const point3d &p0, const point3d &p1,
                        const point3d &p2, const point3d &p3)
{
    return std::max(std::max(dist2(p0,p1), std::max(dist2(p0,p2), dist2(p0,p3))),
                    std::max(dist2(p1,p2), std::max(dist2(p1,p3), dist2(p2,p3))));
}

////////////////////////////////////////////////////////////////////////////////

}  // namespace math
}  // namespace steps

// END

